#!/usr/bin/env python3
"""
validate_all_v2: unified validation of seq_sim_v2.
  Phase 1 (core): topology + branch-length + kappa  — 1215 runs
  Phase 2 (gamma): rate-heterogeneity parameter recovery — 480 runs (FULL)

Usage:  python3 validate_all_v2.py
  SMALL_CORE=True  → 1 tree in core matrix
  SMALL_GAMMA=True → 1 tree in gamma matrix
"""

import sys
import os
import re
import csv
import glob
import time
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

# ============================================================
# path / config
# ============================================================
TREE_DIR = "/mnt/storage5/yfdai/homework/0804/op_stimu_sour/tree"
OUT_DIR  = "/mnt/storage5/yfdai/homework/0805/seq_sim_v2/all_test"
SEQ_SIM_ROOT = "/mnt/storage5/yfdai/homework/0805/"
IQTREE   = "/mnt/storage3/yfdai/download/iqtree-3.1.3-Linux/bin/iqtree3"

# ---- core matrix ----
MODELS       = ["JC69", "HKY", "GTR"]
TREE_HEIGHTS = [0.1, 0.5, 1.0]
SEED_BASE    = 42
NUM_REPLICATES = 3
TIMEOUT_SIM = 300
TIMEOUT_IQ  = 600

IQTREE_MODEL_MAP = {"JC69": "JC", "HKY": "HKY", "GTR": "GTR"}

LENGTHS_BY_SIZE = {
    "small": [500, 1000, 2000],
    "large": [2000, 5000, 10000],
}

MODEL_EXTRA = {
    "JC69": [],
    "HKY":  ["-t", "1.0"],
    "GTR":  ["-r", "1.0", "2.0", "1.0", "1.0", "2.0", "1.0"],
}

# -t 1.0 with equal freqs (0.25): kappa = (1.0 * 0.25) / 0.125 = 2.0
KAPPA_TRUE = 2.0

SMALL_CORE = False

# ---- gamma matrix ----
GAMMA_HEIGHT = 0.5
GAMMA_LENGTHS = [500, 1000, 2000, 5000]

RATE_CONFIGS = [
    ("plain",          None, None, None, "HKY"),
    ("G4_a0.5",        4,    0.5,  None, "HKY+G4"),
    ("G4_a1.0",        4,    1.0,  None, "HKY+G4"),
    ("G4_a2.0",        4,    2.0,  None, "HKY+G4"),
    ("G8_a1.0",        8,    1.0,  None, "HKY+G8"),
    ("I0.2",           None, None, 0.2,  "HKY+I"),
    ("I0.2_G4_a1.0",   4,    1.0,  0.2,  "HKY+I+G4"),
    ("I0.2_G8_a1.0",   8,    1.0,  0.2,  "HKY+I+G8"),
]

SMALL_GAMMA = True

NUM_WORKERS = 16
IQTREE_THREADS = "6"


# ============================================================
# Newick splited / RF
# ============================================================
_NEWICK_RE = re.compile(r':[0-9]+(\.[0-9]+)?([eE][+\-]?[0-9]+)?')


def _strip_branch_lengths(text):
    return _NEWICK_RE.sub('', text)


def _parse_newick_splits(newick_text):
    s = newick_text.strip()
    if s.endswith(';'):
        s = s[:-1]
    s = _strip_branch_lengths(s)
    s = re.sub(r'\s+', '', s)

    all_tips = set()
    raw_splits = []

    def _parse(i):
        nonlocal all_tips, raw_splits
        if i >= len(s):
            return i, set()
        if s[i] == '(':
            i += 1
            leaves = set()
            while True:
                i, child = _parse(i)
                leaves.update(child)
                if i >= len(s):
                    break
                if s[i] == ')':
                    break
                if s[i] == ',':
                    i += 1
            i += 1
            if len(leaves) > 1:
                raw_splits.append(frozenset(leaves))
            return i, leaves
        else:
            j = i
            while j < len(s) and s[j] not in ',()':
                j += 1
            name = s[i:j]
            if name:
                all_tips.add(name)
            return j, {name} if name else set()

    _parse(0)
    full = frozenset(all_tips)
    ref = sorted(all_tips)[0] if all_tips else None

    splits = set()
    for sp in raw_splits:
        if sp == full:
            continue
        complement = frozenset(full - sp)
        if len(sp) <= 1 or len(complement) <= 1:
            continue
        canonical = complement if ref in sp else sp
        splits.add(canonical)
    return full, splits


def compute_rf(newick_a, newick_b):
    tips_a, splits_a = _parse_newick_splits(newick_a)
    tips_b, splits_b = _parse_newick_splits(newick_b)
    if tips_a != tips_b:
        return None, False
    rf = len(splits_a - splits_b) + len(splits_b - splits_a)
    return rf, (rf == 0)


# ============================================================
# Newick branch-length parser
# ============================================================

def _read_float(s, i):
    j = i
    while j < len(s) and (s[j].isdigit() or s[j] in '+-.eE'):
        j += 1
    if j > i:
        return j, float(s[i:j])
    raise ValueError(f"Expected branch length at position {i}, got '{s[i:i+30]}'")


def _parse_tree(s, i):
    if s[i] == '(':
        i += 1
        children = []
        while True:
            i, child = _parse_tree(s, i)
            children.append(child)
            if i >= len(s) or s[i] == ')':
                break
            if s[i] == ',':
                i += 1
        i += 1

        j = i
        while j < len(s) and s[j] not in ':,)();':
            j += 1
        i = j

        length = 0.0
        if i < len(s) and s[i] == ':':
            i, length = _read_float(s, i + 1)

        tips = set()
        for child in children:
            tips.update(child['tips'])

        return i, {'tips': frozenset(tips), 'length': length, 'children': children}
    else:
        j = i
        while j < len(s) and s[j] not in ':,)();':
            j += 1
        name = s[i:j].strip()
        i = j

        length = 0.0
        if i < len(s) and s[i] == ':':
            i, length = _read_float(s, i + 1)

        tips = {name} if name else set()
        return i, {'tips': frozenset(tips), 'length': length, 'children': []}


def _compute_tree_height(root_node):
    max_dist = 0.0

    def _dfs(node, dist):
        nonlocal max_dist
        if not node['children']:
            max_dist = max(max_dist, dist)
        else:
            for child in node['children']:
                _dfs(child, dist + child['length'])

    _dfs(root_node, 0.0)
    return max_dist


def _canonical_edge_key(tips, all_tips, ref_tip):
    side = frozenset(tips)
    complement = frozenset(all_tips - side)
    if not side or not complement:
        return None, None

    if len(side) == 1:
        return ("TIP", next(iter(side))), True
    if len(complement) == 1:
        return ("TIP", next(iter(complement))), True

    canonical = complement if ref_tip in side else side
    return ("SPLIT", tuple(sorted(canonical))), False


def _parse_newick_edges(nwk, target_height=None):
    nwk = nwk.strip()
    if not nwk:
        return {}, frozenset()

    s = nwk
    if s.endswith(';'):
        s = s[:-1]
    s = re.sub(r'\s+', '', s)
    if not s:
        return {}, frozenset()

    pos, root = _parse_tree(s, 0)
    if pos != len(s):
        raise ValueError(
            f"Newick only partially parsed: pos={pos}/{len(s)}, "
            f"remaining: '{s[pos:pos+40]}'"
        )
    all_tips = root['tips']
    if not all_tips:
        raise ValueError("Tree contains no taxa")

    if target_height is not None and target_height > 0:
        tree_height = _compute_tree_height(root)
        if tree_height <= 1e-12:
            raise ValueError("Cannot scale a zero-height tree")
        scale = target_height / tree_height
    else:
        scale = 1.0

    ref_tip = min(all_tips) if all_tips else None
    edges = {}

    def _collect(node):
        for child in node['children']:
            _collect(child)
            tips = child['tips']
            length = child['length'] * scale

            if len(tips) == 0 or len(tips) == len(all_tips):
                continue

            key, is_terminal = _canonical_edge_key(tips, all_tips, ref_tip)
            if key is None:
                continue

            if key in edges:
                edges[key]['length'] += length
            else:
                edges[key] = {
                    'length': length,
                    'is_terminal': is_terminal,
                }

    _collect(root)
    return edges, all_tips


def _compute_branch_errors(true_nwk, inferred_nwk, target_height):
    true_edges, _ = _parse_newick_edges(true_nwk, target_height)
    inf_edges, _ = _parse_newick_edges(inferred_nwk)

    pairs = []
    for key, te in true_edges.items():
        if key in inf_edges:
            pairs.append((te['length'], inf_edges[key]['length']))

    n_true = len(true_edges)
    matched_fraction = len(pairs) / max(n_true, 1)

    if not pairs:
        return {
            'n_matched': 0, 'n_true_total': n_true,
            'matched_fraction': matched_fraction,
            'mae': None, 'rmse': None, 'pearson_r': None, 'mean_rel_error': None,
            'slope': None, 'total_length_ratio': None,
        }

    true_arr = np.array([p[0] for p in pairs], dtype=np.float64)
    infer_arr = np.array([p[1] for p in pairs], dtype=np.float64)

    finite_mask = np.isfinite(true_arr) & np.isfinite(infer_arr)
    if not finite_mask.all():
        true_arr = true_arr[finite_mask]
        infer_arr = infer_arr[finite_mask]

    if len(true_arr) == 0:
        return {
            'n_matched': len(pairs), 'n_true_total': n_true,
            'matched_fraction': matched_fraction,
            'mae': None, 'rmse': None, 'pearson_r': None, 'mean_rel_error': None,
            'slope': None, 'total_length_ratio': None,
        }

    mae = float(np.mean(np.abs(true_arr - infer_arr)))
    rmse = float(np.sqrt(np.mean((true_arr - infer_arr) ** 2)))

    denom_slope = float(np.dot(true_arr, true_arr))
    slope = float(np.dot(true_arr, infer_arr) / denom_slope) if denom_slope > 1e-15 else None
    total_length_ratio = float(infer_arr.sum() / max(true_arr.sum(), 1e-15))

    n = len(true_arr)
    if n >= 3:
        xm = true_arr - true_arr.mean()
        ym = infer_arr - infer_arr.mean()
        denom = np.sqrt(np.sum(xm * xm) * np.sum(ym * ym))
        pearson_r = float(np.sum(xm * ym) / denom) if denom > 1e-15 else None
    else:
        pearson_r = None

    valid_mask = true_arr > 1e-6
    if valid_mask.sum() > 0:
        rel_errs = np.abs(true_arr[valid_mask] - infer_arr[valid_mask]) / true_arr[valid_mask]
        mean_rel_error = float(np.mean(rel_errs) * 100.0)
    else:
        mean_rel_error = None

    return {
        'n_matched': len(pairs),
        'n_true_total': n_true,
        'matched_fraction': matched_fraction,
        'mae': mae,
        'rmse': rmse,
        'pearson_r': pearson_r,
        'mean_rel_error': mean_rel_error,
        'slope': slope,
        'total_length_ratio': total_length_ratio,
    }


# ============================================================
# IQ-TREE .iqtree parameter parser
# ============================================================

def _parse_iqtree_params(iqtree_path):
    result = {"kappa": None, "alpha": None, "p_invar": None, "freq_max_dev": None}
    try:
        with open(iqtree_path, encoding="utf-8", errors="replace") as f:
            text = f.read()
    except (OSError, UnicodeDecodeError):
        return result
    if not text:
        return result

    m = re.search(
        r'(?:Ti/tv ratio|kappa|ts/tv)\s*(?:\([^)]*\))?\s*[:=]\s*([\d.]+)',
        text, re.IGNORECASE,
    )
    if m:
        try:
            result["kappa"] = float(m.group(1))
        except ValueError:
            pass

    if result["kappa"] is None:
        m = re.search(r'A-G:\s*([\d.]+)', text)
        if m:
            try:
                result["kappa"] = float(m.group(1))
            except ValueError:
                pass

    m = re.search(
        r'(?:Gamma shape(?:[\s-]*alpha)?)\s*(?:\([^)]*\))?\s*[:=]\s*([\d.]+)',
        text, re.IGNORECASE,
    )
    if m:
        try:
            result["alpha"] = float(m.group(1))
        except ValueError:
            pass

    m = re.search(
        r'(?:Proportion of invariable sites|p[._-]?invar)\s*(?:\([^)]*\))?\s*[:=]\s*([\d.]+)',
        text, re.IGNORECASE,
    )
    if m:
        try:
            result["p_invar"] = float(m.group(1))
        except ValueError:
            pass

    freqs = {}
    for m in re.finditer(r'pi\(([ACGT])\)\s*=\s*([\d.]+)', text):
        freqs[m.group(1)] = float(m.group(2))
    if len(freqs) == 4:
        result["freq_max_dev"] = max(abs(freqs[b] - 0.25) for b in "ACGT")
    elif len(freqs) > 0:
        devs = [abs(freqs.get(b, 0.25) - 0.25) for b in "ACGT"]
        result["freq_max_dev"] = max(devs)

    return result


# ============================================================
# general helpers
# ============================================================

def _preflight():
    """Verify external dependencies exist before launching tests."""
    checks = [
        (os.path.isdir(TREE_DIR), f"TREE_DIR not found: {TREE_DIR}"),
        (os.path.isfile(IQTREE),  f"IQ-TREE not found: {IQTREE}"),
        (os.access(IQTREE, os.X_OK), f"IQ-TREE not executable: {IQTREE}"),
        (os.path.isdir(SEQ_SIM_ROOT), f"SEQ_SIM_ROOT not found: {SEQ_SIM_ROOT}"),
    ]
    for ok, msg in checks:
        if not ok:
            raise RuntimeError(msg)

    try:
        proc = subprocess.run([IQTREE, "--version"], capture_output=True, text=True, timeout=10)
        return proc.stdout.strip().split("\n")[0] if proc.stdout else "unknown"
    except Exception:
        return "unknown"


def _cleanup_prefix(prefix):
    """Remove existing IQ-TREE / seq_sim output files for a given prefix."""
    for path in glob.glob(prefix + ".*"):
        try:
            os.remove(path)
        except OSError:
            pass


def _make_env():
    """Build subprocess environment with PYTHONPATH and BLAS thread limits."""
    env = os.environ.copy()
    env["PYTHONPATH"] = SEQ_SIM_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    env["OMP_NUM_THREADS"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "1"
    env["MKL_NUM_THREADS"] = "1"
    env["NUMEXPR_NUM_THREADS"] = "1"
    return env


def _split_newick_content(content):
    trees = []
    depth = 0
    start = 0
    for i, ch in enumerate(content):
        if ch == '(':
            depth += 1
        elif ch == ')':
            depth -= 1
        elif ch == ';' and depth == 0:
            tree = content[start:i + 1].strip()
            if tree and '(' in tree:
                trees.append(tree)
            start = i + 1
    remaining = content[start:].strip()
    if remaining and '(' in remaining:
        trees.append(remaining)
    return trees


def _get_lengths_for_tree(basename):
    if "20tip" in basename.lower():
        return LENGTHS_BY_SIZE["large"]
    return LENGTHS_BY_SIZE["small"]


def _dump_error(out_dir, tag, stage, stderr):
    if not stderr:
        return
    path = os.path.join(out_dir, tag, f"{tag}_{stage}_error.log")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(stderr)


# ============================================================
# Phase 1 — core test (topology + blen + kappa)
# ============================================================

def run_core_test(newick_str, model, length, height, tag, seed):
    result = {
        "tag": tag, "model": model, "length": length, "height": height,
        "seed": seed, "status": "OK",
        "rf_distance": None, "topology_match": None,
        "n_matched": None, "n_true_total": None, "matched_fraction": None,
        "mae": None, "rmse": None, "pearson_r": None, "mean_rel_error": None,
        "slope": None, "total_length_ratio": None,
        "rf0_mae": None, "rf0_rmse": None, "rf0_pearson_r": None,
        "rf0_slope": None, "rf0_total_length_ratio": None,
        "kappa_est": None, "kappa_err": None,
        "freq_max_dev": None,
        "time_sim": None, "time_iqtree": None, "error": None,
    }

    os.makedirs(OUT_DIR, exist_ok=True)

    tmp_tree = os.path.join(OUT_DIR, tag, "_tree.nwk")
    os.makedirs(os.path.dirname(tmp_tree), exist_ok=True)
    with open(tmp_tree, "w") as f:
        f.write(newick_str)

    env = _make_env()

    out_phy = os.path.join(OUT_DIR, tag, f"{tag}.phy")
    if os.path.isfile(out_phy):
        os.remove(out_phy)

    # --- seq_sim_v2 ---
    t0 = time.time()
    cmd_sim = [
        sys.executable, "-m", "seq_sim_v2",
        "-m", model,
        "-l", str(length),
        "-d", str(height),
        "-z", str(seed),
        "-o", "p",
        "-y", tag,
        "-q",
        tmp_tree,
        *MODEL_EXTRA.get(model, []),
    ]
    try:
        proc = subprocess.run(
            cmd_sim, cwd=OUT_DIR, env=env,
            capture_output=True, text=True, timeout=TIMEOUT_SIM,
        )
    except subprocess.TimeoutExpired:
        result["status"] = "FAIL"; result["error"] = "seq_sim_v2 timeout"
        return result
    result["time_sim"] = round(time.time() - t0, 1)

    if proc.returncode != 0:
        result["status"] = "FAIL"
        _dump_error(OUT_DIR, tag, "seq_sim_v2", proc.stderr)
        result["error"] = f"seq_sim_v2 rc={proc.returncode}: {proc.stderr[:200]}"
        return result

    if not os.path.isfile(out_phy):
        result["status"] = "FAIL"
        result["error"] = f"seq_sim_v2 output not found: {out_phy}"
        return result

    # --- IQ-TREE ---
    iq_model = IQTREE_MODEL_MAP.get(model, model)
    iq_prefix = os.path.join(OUT_DIR, tag, f"{tag}_iqtree")
    _cleanup_prefix(iq_prefix)
    t0 = time.time()
    cmd_iq = [
        IQTREE, "-s", out_phy, "-m", iq_model,
        "--prefix", iq_prefix, "-T", IQTREE_THREADS, "--quiet",
    ]
    try:
        proc = subprocess.run(cmd_iq, capture_output=True, text=True,
                              timeout=TIMEOUT_IQ, env=env)
    except subprocess.TimeoutExpired:
        result["status"] = "FAIL"; result["error"] = "IQ-TREE timeout"
        return result
    result["time_iqtree"] = round(time.time() - t0, 1)

    if proc.returncode != 0:
        result["status"] = "FAIL"
        _dump_error(OUT_DIR, tag, "iqtree", proc.stderr)
        result["error"] = f"IQ-TREE rc={proc.returncode}: {proc.stderr[:200]}"
        return result

    iq_treefile = os.path.join(OUT_DIR, tag, f"{tag}_iqtree.treefile")
    if not os.path.isfile(iq_treefile):
        result["status"] = "FAIL"
        result["error"] = f"IQ-TREE treefile not found: {iq_treefile}"
        return result

    with open(iq_treefile) as f:
        inferred = f.read().strip()

    # --- RF ---
    rf_dist, match = compute_rf(newick_str, inferred)
    result["rf_distance"] = rf_dist
    result["topology_match"] = match

    # --- Branch length ---
    blen = _compute_branch_errors(newick_str, inferred, height)
    result["n_matched"] = blen["n_matched"]
    result["n_true_total"] = blen["n_true_total"]
    result["matched_fraction"] = blen["matched_fraction"]
    result["mae"] = blen["mae"]
    result["rmse"] = blen["rmse"]
    result["pearson_r"] = blen["pearson_r"]
    result["mean_rel_error"] = blen["mean_rel_error"]
    result["slope"] = blen["slope"]
    result["total_length_ratio"] = blen["total_length_ratio"]

    if result["topology_match"]:
        result["rf0_mae"] = blen["mae"]
        result["rf0_rmse"] = blen["rmse"]
        result["rf0_pearson_r"] = blen["pearson_r"]
        result["rf0_slope"] = blen["slope"]
        result["rf0_total_length_ratio"] = blen["total_length_ratio"]

    # --- kappa + freqs (HKY only) ---
    if model == "HKY":
        iqtree_file = os.path.join(OUT_DIR, tag, f"{tag}_iqtree.iqtree")
        params = _parse_iqtree_params(iqtree_file)
        result["kappa_est"] = params["kappa"]
        result["freq_max_dev"] = params["freq_max_dev"]
        if params["kappa"] is not None:
            result["kappa_err"] = abs(params["kappa"] - KAPPA_TRUE)
        else:
            result["status"] = "FAIL"
            result["error"] = "Failed to parse HKY kappa from .iqtree"

    return result


# ============================================================
# Phase 1 aggregation
# ============================================================

def _aggregate_core(results):
    groups = defaultdict(list)
    for r in results:
        base = re.sub(r'_r\d+$', '', r["tag"])
        key = (base, r["model"], r["length"], r["height"])
        groups[key].append(r)

    topo_agg, blen_agg = [], []
    for (base, model, length, height), reps in groups.items():
        ok = [r for r in reps if r["status"] == "OK"]
        n_ok = len(ok)
        n_expected = NUM_REPLICATES
        complete = n_ok == n_expected

        rf_vals = [r["rf_distance"] for r in ok if r.get("rf_distance") is not None]
        pass_count = sum(1 for r in ok if r.get("topology_match"))
        majority = complete and pass_count >= (n_expected // 2 + 1)

        topo_agg.append({
            "tag": base, "model": model, "length": length, "height": height,
            "replicates_expected": n_expected,
            "replicates_ok": n_ok,
            "complete": complete,
            "avg_rf": round(float(np.mean(rf_vals)), 1) if rf_vals else None,
            "min_rf": min(rf_vals) if rf_vals else None,
            "max_rf": max(rf_vals) if rf_vals else None,
            "rf_pass_count": pass_count,
            "majority_topo_match": majority,
            "failed": len(reps) - n_ok,
        })

        mae_v  = [r["mae"] for r in ok if r.get("mae") is not None]
        rmse_v = [r["rmse"] for r in ok if r.get("rmse") is not None]
        pr_v   = [r["pearson_r"] for r in ok if r.get("pearson_r") is not None]
        rel_v  = [r["mean_rel_error"] for r in ok if r.get("mean_rel_error") is not None]
        mt_v   = [r["n_matched"] for r in ok if r.get("n_matched") is not None]
        mf_v   = [r["matched_fraction"] for r in ok if r.get("matched_fraction") is not None]
        sl_v   = [r["slope"] for r in ok if r.get("slope") is not None]
        tlr_v  = [r["total_length_ratio"] for r in ok if r.get("total_length_ratio") is not None]
        ke_v   = [r["kappa_err"] for r in ok if r.get("kappa_err") is not None]
        fd_v   = [r["freq_max_dev"] for r in ok if r.get("freq_max_dev") is not None]

        r0m_v  = [r["rf0_mae"] for r in ok if r.get("rf0_mae") is not None]
        r0r_v  = [r["rf0_rmse"] for r in ok if r.get("rf0_rmse") is not None]
        r0p_v  = [r["rf0_pearson_r"] for r in ok if r.get("rf0_pearson_r") is not None]
        r0s_v  = [r["rf0_slope"] for r in ok if r.get("rf0_slope") is not None]
        r0t_v  = [r["rf0_total_length_ratio"] for r in ok if r.get("rf0_total_length_ratio") is not None]

        def _fisher_mean(vals):
            if not vals:
                return None
            x = np.clip(np.asarray(vals, dtype=np.float64), -0.999999, 0.999999)
            return round(float(np.tanh(np.mean(np.arctanh(x)))), 4)

        fisher_mean_r = _fisher_mean(pr_v)
        rf0_fisher_mean_r = _fisher_mean(r0p_v)

        blen_agg.append({
            "tag": base, "model": model, "length": length, "height": height,
            "replicates_expected": n_expected,
            "replicates_ok": n_ok,
            "complete": complete,
            "mean_mae": round(float(np.mean(mae_v)), 6) if mae_v else None,
            "mean_rmse": round(float(np.mean(rmse_v)), 6) if rmse_v else None,
            "mean_pearson_r": fisher_mean_r,
            "mean_rel_error_pct": round(float(np.mean(rel_v)), 2) if rel_v else None,
            "mean_n_matched": round(float(np.mean(mt_v)), 1) if mt_v else None,
            "mean_matched_fraction": round(float(np.mean(mf_v)), 4) if mf_v else None,
            "mean_slope": round(float(np.mean(sl_v)), 4) if sl_v else None,
            "mean_total_length_ratio": round(float(np.mean(tlr_v)), 4) if tlr_v else None,
            "mean_rf0_mae": round(float(np.mean(r0m_v)), 6) if r0m_v else None,
            "mean_rf0_rmse": round(float(np.mean(r0r_v)), 6) if r0r_v else None,
            "mean_rf0_pearson_r": rf0_fisher_mean_r,
            "mean_rf0_slope": round(float(np.mean(r0s_v)), 4) if r0s_v else None,
            "mean_rf0_total_length_ratio": round(float(np.mean(r0t_v)), 4) if r0t_v else None,
            "rf0_count": len(r0m_v),
            "mean_kappa_err": round(float(np.mean(ke_v)), 6) if ke_v else None,
            "mean_freq_max_dev": round(float(np.mean(fd_v)), 6) if fd_v else None,
            "rf_pass_count": pass_count,
            "majority_topo_match": majority,
            "failed": len(reps) - n_ok,
        })

    return topo_agg, blen_agg


# ============================================================
# Phase 2 — gamma test
# ============================================================

def run_gamma_test(newick_str, cfg_label, g_cats, alpha, p_invar,
                   iq_model, length, tag, seed):
    result = {
        "tag": tag, "config_label": cfg_label, "length": length,
        "seed": seed, "status": "OK",
        "alpha_true": alpha, "kappa_true": KAPPA_TRUE, "p_invar_true": p_invar,
        "alpha_est": None, "kappa_est": None, "p_invar_est": None,
        "alpha_err": None, "kappa_err": None, "p_invar_err": None,
        "freq_max_dev": None,
        "time_sim": None, "time_iqtree": None, "error": None,
    }

    os.makedirs(OUT_DIR, exist_ok=True)

    tmp_tree = os.path.join(OUT_DIR, tag, "_tree.nwk")
    os.makedirs(os.path.dirname(tmp_tree), exist_ok=True)
    with open(tmp_tree, "w") as f:
        f.write(newick_str)

    env = _make_env()

    out_phy = os.path.join(OUT_DIR, tag, f"{tag}.phy")
    if os.path.isfile(out_phy):
        os.remove(out_phy)

    # --- seq_sim_v2 ---
    extra_sim = []
    if g_cats is not None:
        extra_sim.extend(["-g", str(g_cats)])
    if alpha is not None:
        extra_sim.extend(["-a", str(alpha)])
    if p_invar is not None:
        extra_sim.extend(["-i", str(p_invar)])

    cmd_sim = [
        sys.executable, "-m", "seq_sim_v2",
        "-m", "HKY",
        "-l", str(length),
        "-d", str(GAMMA_HEIGHT),
        "-z", str(seed),
        "-t", "1.0",
        "-o", "p",
        "-y", tag,
        "-q",
        tmp_tree,
        *extra_sim,
    ]

    t0 = time.time()
    try:
        proc = subprocess.run(cmd_sim, cwd=OUT_DIR, env=env,
                              capture_output=True, text=True, timeout=TIMEOUT_SIM)
    except subprocess.TimeoutExpired:
        result["status"] = "FAIL"; result["error"] = "seq_sim_v2 timeout"
        return result
    result["time_sim"] = round(time.time() - t0, 1)

    if proc.returncode != 0:
        result["status"] = "FAIL"
        _dump_error(OUT_DIR, tag, "seq_sim_v2", proc.stderr)
        result["error"] = f"seq_sim_v2 rc={proc.returncode}: {proc.stderr[:200]}"
        return result

    if not os.path.isfile(out_phy):
        result["status"] = "FAIL"
        result["error"] = f"seq_sim_v2 output not found: {out_phy}"
        return result

    # --- IQ-TREE ---
    iq_prefix = os.path.join(OUT_DIR, tag, f"{tag}_iqtree")
    _cleanup_prefix(iq_prefix)
    t0 = time.time()
    cmd_iq = [
        IQTREE, "-s", out_phy, "-m", iq_model,
        "--prefix", iq_prefix, "-T", IQTREE_THREADS, "--quiet",
    ]
    try:
        proc = subprocess.run(cmd_iq, capture_output=True, text=True,
                              timeout=TIMEOUT_IQ, env=env)
    except subprocess.TimeoutExpired:
        result["status"] = "FAIL"; result["error"] = "IQ-TREE timeout"
        return result
    result["time_iqtree"] = round(time.time() - t0, 1)

    if proc.returncode != 0:
        result["status"] = "FAIL"
        _dump_error(OUT_DIR, tag, "iqtree", proc.stderr)
        result["error"] = f"IQ-TREE rc={proc.returncode}: {proc.stderr[:200]}"
        return result

    iqtree_file = os.path.join(OUT_DIR, tag, f"{tag}_iqtree.iqtree")
    if not os.path.isfile(iqtree_file):
        result["status"] = "FAIL"
        result["error"] = f"IQ-TREE .iqtree not found: {iqtree_file}"
        return result

    params = _parse_iqtree_params(iqtree_file)
    result["kappa_est"] = params["kappa"]
    result["alpha_est"] = params["alpha"]
    result["p_invar_est"] = params["p_invar"]
    result["freq_max_dev"] = params["freq_max_dev"]

    missing = []
    if params["kappa"] is None:
        missing.append("kappa")
    if alpha is not None and params["alpha"] is None:
        missing.append("alpha")
    if p_invar is not None and params["p_invar"] is None:
        missing.append("p_invar")
    if missing:
        result["status"] = "FAIL"
        result["error"] = "IQ-TREE parameter parsing failed: " + ", ".join(missing)
        return result

    bad = []
    if result["p_invar_est"] is not None and (result["p_invar_est"] > 1.0
                                              or result["p_invar_est"] < 0.0):
        bad.append(f"p_invar={result['p_invar_est']:.4f}")
    if result["alpha_est"] is not None and result["alpha_est"] > 50.0:
        bad.append(f"alpha={result['alpha_est']:.1f}")
    if bad:
        result["status"] = "FAIL"
        result["error"] = "IQ-TREE estimate out of physical range: " + ", ".join(bad)
        return result

    if result["kappa_est"] is not None:
        result["kappa_err"] = abs(result["kappa_est"] - KAPPA_TRUE)
    if alpha is not None and result["alpha_est"] is not None:
        result["alpha_err"] = abs(result["alpha_est"] - alpha)
    if p_invar is not None and result["p_invar_est"] is not None:
        result["p_invar_err"] = abs(result["p_invar_est"] - p_invar)

    return result


def _aggregate_gamma(results):
    groups = defaultdict(list)
    for r in results:
        base = re.sub(r'_r\d+$', '', r["tag"])
        key = (base, r["config_label"], r["length"])
        groups[key].append(r)

    aggregated = []
    for (base, cfg, length), reps in groups.items():
        ok = [r for r in reps if r["status"] == "OK"]
        n_ok = len(ok)
        n_expected = NUM_REPLICATES
        complete = n_ok == n_expected

        ae = [r["alpha_err"] for r in ok if r.get("alpha_err") is not None]
        ke = [r["kappa_err"] for r in ok if r.get("kappa_err") is not None]
        pe = [r["p_invar_err"] for r in ok if r.get("p_invar_err") is not None]
        fd = [r["freq_max_dev"] for r in ok if r.get("freq_max_dev") is not None]

        a_est = [r["alpha_est"] for r in ok if r.get("alpha_est") is not None]
        k_est = [r["kappa_est"] for r in ok if r.get("kappa_est") is not None]
        p_est = [r["p_invar_est"] for r in ok if r.get("p_invar_est") is not None]

        alpha_ref = next((r.get("alpha_true") for r in reps if r.get("alpha_true") is not None), None)
        kappa_ref = next((r.get("kappa_true") for r in reps if r.get("kappa_true") is not None), None)
        pinv_ref = next((r.get("p_invar_true") for r in reps if r.get("p_invar_true") is not None), None)

        def _bias_sd(ests, true_val):
            if not ests or true_val is None:
                return None, None
            x = np.asarray(ests, dtype=np.float64)
            bias = float(np.mean(x)) - true_val
            sd = float(np.std(x, ddof=1)) if len(x) >= 2 else None
            return round(bias, 6), round(sd, 6) if sd is not None else None

        a_bias, a_sd = _bias_sd(a_est, alpha_ref)
        k_bias, k_sd = _bias_sd(k_est, kappa_ref)
        p_bias, p_sd = _bias_sd(p_est, pinv_ref)

        aggregated.append({
            "tag": base,
            "config_label": cfg,
            "length": length,
            "replicates_expected": n_expected,
            "replicates_ok": n_ok,
            "complete": complete,
            "mean_alpha_est": round(float(np.mean(a_est)), 4) if a_est else None,
            "mean_alpha_err": round(float(np.mean(ae)), 6) if ae else None,
            "alpha_bias": a_bias,
            "alpha_sd": a_sd,
            "mean_kappa_est": round(float(np.mean(k_est)), 4) if k_est else None,
            "mean_kappa_err": round(float(np.mean(ke)), 6) if ke else None,
            "kappa_bias": k_bias,
            "kappa_sd": k_sd,
            "mean_p_invar_est": round(float(np.mean(p_est)), 4) if p_est else None,
            "mean_p_invar_err": round(float(np.mean(pe)), 6) if pe else None,
            "p_invar_bias": p_bias,
            "p_invar_sd": p_sd,
            "mean_freq_max_dev": round(float(np.mean(fd)), 6) if fd else None,
            "failed": len(reps) - n_ok,
        })
    return aggregated


# ============================================================
# main
# ============================================================

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    iqtree_ver = _preflight()
    print(f"IQ-TREE: {iqtree_ver}")
    print()

    # --- Discover tree files ---
    tree_files = sorted(glob.glob(os.path.join(TREE_DIR, "*.nwk")))
    if not tree_files:
        tree_files = sorted(glob.glob(os.path.join(TREE_DIR, "*.tree")))
    if not tree_files:
        tree_files = sorted(glob.glob(os.path.join(TREE_DIR, "*.tre")))

    print(f"== Found {len(tree_files)} tree file(s)")
    for tf in tree_files:
        print(f"   {os.path.basename(tf)}")

    if not tree_files:
        raise RuntimeError(f"No .nwk/.tree/.tre files found in {TREE_DIR}")

    all_trees = []
    for tf in tree_files:
        with open(tf) as f:
            content = f.read()
        individual = _split_newick_content(content)
        for idx, nwk in enumerate(individual):
            all_trees.append((os.path.basename(tf), nwk, idx, len(individual)))
    print(f"   total individual trees: {len(all_trees)}")
    if not all_trees:
        raise RuntimeError("No valid semicolon-terminated Newick trees parsed from tree files")

    # ================================================================
    # quick parser sanity check
    t1 = "((A:1,B:1):0.2,(C:1,D:1):0.3);"
    t2 = "(A:1,B:1,(C:1,D:1):0.5);"
    rf_test, match_test = compute_rf(t1, t2)
    e1, tips1 = _parse_newick_edges(t1)
    e2, tips2 = _parse_newick_edges(t2)
    key_int = ("SPLIT", ("C", "D"))
    if rf_test != 0 or not match_test:
        raise RuntimeError(f"Self-test RF failed: rf={rf_test}, match={match_test}")
    if tips1 != tips2:
        raise RuntimeError("Self-test taxon sets differ")
    if key_int not in e1 or not np.isclose(e1[key_int]["length"], 0.5, atol=1e-12):
        raise RuntimeError(f"Self-test e1 internal edge: {e1.get(key_int)}")
    if key_int not in e2 or not np.isclose(e2[key_int]["length"], 0.5, atol=1e-12):
        raise RuntimeError(f"Self-test e2 internal edge: {e2.get(key_int)}")

    t3 = "(A:0.2,(B:1,C:1,D:1):0.3);"
    t4 = "(A:0.5,B:1,C:1,D:1);"
    e3, tips3 = _parse_newick_edges(t3)
    e4, tips4 = _parse_newick_edges(t4)
    key_tip = ("TIP", "A")
    if tips3 != tips4:
        raise RuntimeError("Self-test singleton-root taxon sets differ")
    if key_tip not in e3 or not np.isclose(e3[key_tip]["length"], 0.5, atol=1e-12):
        raise RuntimeError(f"Self-test root-adjacent terminal edge: {e3.get(key_tip)}")
    if key_tip not in e4 or not np.isclose(e4[key_tip]["length"], 0.5, atol=1e-12):
        raise RuntimeError(f"Self-test unrooted terminal edge: {e4.get(key_tip)}")
    print("Parser self-test: PASS\n")
    # ================================================================

    core_trees = list(all_trees)
    if SMALL_CORE:
        four_tip = [(b, n, i, t) for b, n, i, t in core_trees if "4tip" in b.lower()]
        core_trees = four_tip[:1] if four_tip else core_trees[:1]
        print(f"\n  SMALL_CORE: using {core_trees[0][0]} tree #{core_trees[0][2]}")

    core_tests = []
    for bname, nwk, idx, total_in_file in core_trees:
        lengths = _get_lengths_for_tree(bname)
        tname = os.path.splitext(bname)[0]
        base_tag = f"{tname}_t{idx}" if total_in_file > 1 else tname
        for model in MODELS:
            for length in lengths:
                for height in TREE_HEIGHTS:
                    for rep in range(NUM_REPLICATES):
                        sd = SEED_BASE + rep
                        tag = f"{base_tag}__{model}_l{length}_h{height}_r{sd}"
                        core_tests.append((nwk, model, length, height, tag, sd))

    n_core = len(core_tests)
    print(f"\n{'='*60}")
    print(f"  Phase 1 — Core: {n_core} combinations")
    print(f"  Models: {MODELS}   Heights: {TREE_HEIGHTS}")
    print(f"  Lengths 4/6t: {LENGTHS_BY_SIZE['small']}   Lengths 20t: {LENGTHS_BY_SIZE['large']}")
    print(f"  Replicates: {NUM_REPLICATES}   Workers: {NUM_WORKERS}")
    print(f"{'='*60}")

    core_results = []
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as ex:
        fm = {}
        for i, (nwk, model, length, height, tag, sd) in enumerate(core_tests):
            f = ex.submit(run_core_test, nwk, model, length, height, tag, sd)
            fm[f] = {"idx": i, "tag": tag, "model": model,
                     "length": length, "height": height, "seed": sd}

        n_done = 0
        for f in as_completed(fm):
            meta = fm[f]; n_done += 1
            try:
                r = f.result()
            except Exception as exc:
                r = {
                    "tag": meta["tag"], "model": meta["model"],
                    "length": meta["length"], "height": meta["height"],
                    "seed": meta["seed"], "status": "FAIL",
                    "error": f"worker: {type(exc).__name__}: {exc}",
                    "rf_distance": None, "topology_match": None,
                    "n_matched": None, "n_true_total": None, "matched_fraction": None,
                    "mae": None, "rmse": None, "pearson_r": None, "mean_rel_error": None,
                    "slope": None, "total_length_ratio": None,
                    "rf0_mae": None, "rf0_rmse": None, "rf0_pearson_r": None,
                    "rf0_slope": None, "rf0_total_length_ratio": None,
                    "kappa_est": None, "kappa_err": None, "freq_max_dev": None,
                    "time_sim": None, "time_iqtree": None,
                }
            core_results.append((meta["idx"], r))

            ct = meta["tag"]
            if r["status"] == "OK":
                mae_s = f"MAE={r['mae']:.6f}" if r['mae'] is not None else "MAE=?"
                topo_s = "MATCH" if r["topology_match"] else ("N/A" if r["topology_match"] is None else "MISMATCH")
                print(f"[{n_done:3d}/{n_core}] {ct}  OK  {mae_s}  topo={topo_s}  "
                      f"sim={r['time_sim']}s  iq={r['time_iqtree']}s")
            else:
                print(f"[{n_done:3d}/{n_core}] {ct}  FAIL  {r['error']}")

    core_results.sort(key=lambda x: x[0])
    core_results = [r for _, r in core_results]

    # --- Core CSV ---
    raw_core_csv = os.path.join(OUT_DIR, "raw_core.csv")
    raw_core_fields = [
        "tag", "model", "length", "height", "seed", "status",
        "rf_distance", "topology_match",
        "n_matched", "n_true_total", "matched_fraction",
        "mae", "rmse", "pearson_r", "mean_rel_error",
        "slope", "total_length_ratio",
        "rf0_mae", "rf0_rmse", "rf0_pearson_r",
        "rf0_slope", "rf0_total_length_ratio",
        "kappa_est", "kappa_err", "freq_max_dev",
        "time_sim", "time_iqtree", "error",
    ]
    with open(raw_core_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=raw_core_fields)
        w.writeheader()
        for r in core_results:
            w.writerow({k: r.get(k, "") for k in raw_core_fields})

    topo_agg, blen_agg = _aggregate_core(core_results)

    topo_csv = os.path.join(OUT_DIR, "topo.csv")
    topo_fields = [
        "tag", "model", "length", "height",
        "replicates_expected", "replicates_ok", "complete",
        "avg_rf", "min_rf", "max_rf",
        "rf_pass_count", "majority_topo_match", "failed",
    ]
    with open(topo_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=topo_fields)
        w.writeheader()
        for a in topo_agg:
            w.writerow(a)

    blen_csv = os.path.join(OUT_DIR, "blen.csv")
    blen_fields = [
        "tag", "model", "length", "height",
        "replicates_expected", "replicates_ok", "complete",
        "mean_mae", "mean_rmse", "mean_pearson_r",
        "mean_rel_error_pct", "mean_n_matched", "mean_matched_fraction",
        "mean_slope", "mean_total_length_ratio",
        "mean_rf0_mae", "mean_rf0_rmse", "mean_rf0_pearson_r",
        "mean_rf0_slope", "mean_rf0_total_length_ratio", "rf0_count",
        "mean_kappa_err", "mean_freq_max_dev",
        "rf_pass_count", "majority_topo_match", "failed",
    ]
    with open(blen_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=blen_fields)
        w.writeheader()
        for a in blen_agg:
            w.writerow(a)

    # --- Core summary ---
    ok_topo = [a for a in topo_agg if a["complete"]]
    match_topo = [a for a in ok_topo if a["majority_topo_match"]]
    print(f"\n  Phase 1 — Topology: {len(match_topo)}/{len(ok_topo)} majority PASS (complete groups only)")

    ok_blen = [a for a in blen_agg if a["replicates_ok"] > 0]
    m_l = [a["mean_mae"] for a in ok_blen if a.get("mean_mae") is not None]
    p_l = [a["mean_pearson_r"] for a in ok_blen if a.get("mean_pearson_r") is not None]
    sl_l = [a["mean_slope"] for a in ok_blen if a.get("mean_slope") is not None]
    if m_l:
        print(f"  Phase 1 — Branch length: mean MAE={np.mean(m_l):.6f}")
    if p_l:
        print(f"    mean Pearson r={np.mean(p_l):.4f}")
    if sl_l:
        print(f"    mean slope={np.mean(sl_l):.4f}")

    ke_l = [a["mean_kappa_err"] for a in ok_blen if a.get("mean_kappa_err") is not None]
    if ke_l:
        print(f"  Phase 1 — HKY kappa err:  {np.mean(ke_l):.6f}")

    # ================================================================
    # Phase 2 — Gamma matrix
    # ================================================================
    gamma_trees = [(b, n, i, t) for b, n, i, t in all_trees if "4tip" in b.lower()]
    if not gamma_trees:
        raise RuntimeError("Gamma phase requires at least one tree with '4tip' in filename")
    if SMALL_GAMMA:
        gamma_trees = gamma_trees[:1]
        print(f"\n  SMALL_GAMMA: using {gamma_trees[0][0]} tree #{gamma_trees[0][2]}")
    else:
        gamma_trees = gamma_trees[:5]
        print(f"\n  FULL_GAMMA: using {len(gamma_trees)} 4-tip trees")

    gamma_tests = []
    for bname, nwk, idx, total_in_file in gamma_trees:
        tname = os.path.splitext(bname)[0]
        base_tag = f"{tname}_t{idx}" if total_in_file > 1 else tname
        for cfg_label, g_cats, alpha, p_invar, iq_model in RATE_CONFIGS:
            for length in GAMMA_LENGTHS:
                for rep in range(NUM_REPLICATES):
                    sd = SEED_BASE + rep
                    tag = (f"{base_tag}__{cfg_label}"
                           f"_l{length}_h{GAMMA_HEIGHT}_r{sd}")
                    gamma_tests.append((
                        nwk, cfg_label, g_cats, alpha, p_invar,
                        iq_model, length, tag, sd,
                    ))

    n_gamma = len(gamma_tests)
    print(f"\n{'='*60}")
    print(f"  Phase 2 — Gamma: {n_gamma} combinations")
    print(f"  Model: HKY (kappa={KAPPA_TRUE})   Height: {GAMMA_HEIGHT}")
    print(f"  Lengths: {GAMMA_LENGTHS}   Rate configs: {len(RATE_CONFIGS)}")
    for lbl, gc, al, pi, im in RATE_CONFIGS:
        print(f"    {lbl:16s}  ->  IQ-TREE: {im}")
    print(f"  Replicates: {NUM_REPLICATES}   Workers: {NUM_WORKERS}")
    print(f"{'='*60}")

    gamma_results = []
    with ThreadPoolExecutor(max_workers=NUM_WORKERS) as ex:
        fm = {}
        for i, (nwk, cfg_label, g_cats, alpha, p_invar,
                iq_model, length, tag, sd) in enumerate(gamma_tests):
            f = ex.submit(run_gamma_test, nwk, cfg_label, g_cats, alpha,
                          p_invar, iq_model, length, tag, sd)
            fm[f] = {"idx": i, "tag": tag, "cfg_label": cfg_label,
                     "length": length, "seed": sd,
                     "alpha_true": alpha,
                     "kappa_true": KAPPA_TRUE,
                     "p_invar_true": p_invar}

        n_done = 0
        for f in as_completed(fm):
            meta = fm[f]; n_done += 1
            try:
                r = f.result()
            except Exception as exc:
                r = {
                    "tag": meta["tag"], "config_label": meta["cfg_label"],
                    "length": meta["length"], "seed": meta["seed"],
                    "status": "FAIL", "error": f"worker: {type(exc).__name__}: {exc}",
                    "alpha_true": meta["alpha_true"],
                    "kappa_true": meta["kappa_true"],
                    "p_invar_true": meta["p_invar_true"],
                    "alpha_est": None, "kappa_est": None, "p_invar_est": None,
                    "alpha_err": None, "kappa_err": None, "p_invar_err": None,
                    "freq_max_dev": None,
                    "time_sim": None, "time_iqtree": None,
                }
            gamma_results.append((meta["idx"], r))

            gt = meta["tag"]
            if r["status"] == "OK":
                bits = []
                if r["alpha_err"] is not None:
                    bits.append(f"alpha_err={r['alpha_err']:.4f}")
                if r["kappa_err"] is not None:
                    bits.append(f"kappa_err={r['kappa_err']:.4f}")
                if r["p_invar_err"] is not None:
                    bits.append(f"pinv_err={r['p_invar_err']:.4f}")
                bi = "  ".join(bits) if bits else "params_est=N/A"
                print(f"[{n_done:3d}/{n_gamma}] {gt}  OK  {bi}  "
                      f"sim={r['time_sim']}s  iq={r['time_iqtree']}s")
            else:
                print(f"[{n_done:3d}/{n_gamma}] {gt}  FAIL  {r['error']}")

    gamma_results.sort(key=lambda x: x[0])
    gamma_results = [r for _, r in gamma_results]

    # --- Gamma CSV ---
    raw_gamma_csv = os.path.join(OUT_DIR, "raw_gamma.csv")
    raw_gamma_fields = [
        "tag", "config_label", "length", "seed", "status",
        "alpha_true", "alpha_est", "alpha_err",
        "kappa_true", "kappa_est", "kappa_err",
        "p_invar_true", "p_invar_est", "p_invar_err",
        "freq_max_dev",
        "time_sim", "time_iqtree", "error",
    ]
    with open(raw_gamma_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=raw_gamma_fields)
        w.writeheader()
        for r in gamma_results:
            w.writerow({k: r.get(k, "") for k in raw_gamma_fields})

    gamma_agg = _aggregate_gamma(gamma_results)

    gamma_csv = os.path.join(OUT_DIR, "gamma.csv")
    gamma_fields = [
        "tag", "config_label", "length",
        "replicates_expected", "replicates_ok", "complete",
        "mean_alpha_est", "mean_alpha_err", "alpha_bias", "alpha_sd",
        "mean_kappa_est", "mean_kappa_err", "kappa_bias", "kappa_sd",
        "mean_p_invar_est", "mean_p_invar_err", "p_invar_bias", "p_invar_sd",
        "mean_freq_max_dev", "failed",
    ]
    with open(gamma_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=gamma_fields)
        w.writeheader()
        for a in gamma_agg:
            w.writerow(a)

    # --- Gamma summary ---
    ok_g = [a for a in gamma_agg if a["replicates_ok"] > 0]
    print()
    print("=" * 60)
    for pn, pk_est, pk_err, pk_bias in [
        ("alpha", "mean_alpha_est", "mean_alpha_err", "alpha_bias"),
        ("kappa", "mean_kappa_est", "mean_kappa_err", "kappa_bias"),
        ("p_invar", "mean_p_invar_est", "mean_p_invar_err", "p_invar_bias")]:
        errs = [a[pk_err] for a in ok_g if a.get(pk_err) is not None]
        biases = [a[pk_bias] for a in ok_g if a.get(pk_bias) is not None]
        if errs:
            bi_s = f"mean bias={np.mean(biases):.6f}" if biases else ""
            print(f"  mean {pn}_err: {np.mean(errs):.6f} "
                  f"(best={min(errs):.6f} worst={max(errs):.6f})  {bi_s}")

    # Alpha recovery across true values
    print()
    for at in [0.5, 1.0, 2.0]:
        ests = [r["alpha_est"] for r in gamma_results
                if r["status"] == "OK" and r["alpha_true"] == at and r["alpha_est"] is not None]
        if ests:
            print(f"  alpha_true={at}: mean_est={np.mean(ests):.4f}  "
                  f"mean_err={np.mean([abs(e - at) for e in ests]):.4f}  n={len(ests)}")

    # G4 vs G8
    g4 = [r["alpha_est"] for r in gamma_results
          if r["status"] == "OK" and r.get("config_label") == "G4_a1.0"
          and r["alpha_est"] is not None]
    g8 = [r["alpha_est"] for r in gamma_results
          if r["status"] == "OK" and r.get("config_label") == "G8_a1.0"
          and r["alpha_est"] is not None]
    if g4 and g8:
        print(f"\n  G4 vs G8 (alpha=1.0): G4 mean={np.mean(g4):.4f}  G8 mean={np.mean(g8):.4f}")

    # I+G dual recovery
    ig = [r for r in gamma_results
          if r["status"] == "OK"
          and r.get("config_label") in ("I0.2_G4_a1.0", "I0.2_G8_a1.0")
          and r["alpha_est"] is not None and r["p_invar_est"] is not None]
    if ig:
        ae = [abs(r["alpha_est"] - 1.0) for r in ig]
        pe = [abs(r["p_invar_est"] - 0.2) for r in ig]
        print(f"  I+G dual recovery: mean alpha_err={np.mean(ae):.4f}  mean pinv_err={np.mean(pe):.4f}")

    print()
    print(f"  Core outputs:  raw_core.csv, topo.csv, blen.csv")
    print(f"  Gamma outputs: raw_gamma.csv, gamma.csv")
    print("=" * 60)


if __name__ == "__main__":
    main()
