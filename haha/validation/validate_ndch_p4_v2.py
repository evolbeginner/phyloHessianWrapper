#!/usr/bin/env python3
"""
validate_ndch_p4_v2.py

Validates seq_sim_v2 NDRH/NDCH correctness using p4 for per-node parameter recovery.

Pipeline:
  seq_sim_v2 (heterogeneous sim) → .nex alignment
  p4 (same hetero model)         → optLogLike() per-node comp/blen estimation
  compare true vs estimated       → MAE, Pearson r, composition error

Usage (on Linux server with compiled p4):
    python3 validate_ndch_p4_v2.py

Adjust P4_BIN, PROJ_ROOT, OUT_DIR below before running.
"""

import subprocess, os, sys, re, json, shutil
import numpy as np

# ==================== Configuration ====================

P4_BIN = "/mnt/storage3/yfdai/download/p4-phylogenetics-master/bin/p4"
OUT_DIR = "./p4_ndch_validation_v2"

# New project root (directory containing validate.py / __init__.py / seq_sim.py).
# Auto-detect from this script's location.
PROJ_ROOT = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
PROJ_PARENT = os.path.dirname(PROJ_ROOT)

# Test tree (4-tip rooted, clade-level NDCH):
#   AB clade (A,B): f=1 → alt composition
#   CD clade (C,D): f=0 → equal [0.25,0.25,0.25,0.25]
#   Root clade:     f=0 → equal
#   All external branches: 0.2, internal branches: 0.05
CLEAN_TREE = "((A:0.2,B:0.2):0.05,(C:0.2,D:0.2):0.05);"
ANNOTATED_TREE = (
    "((A:0.2[f=1],B:0.2[f=1]):0.05[f=0],"
    "(C:0.2[f=0],D:0.2[f=0]):0.05[f=0]);"
)

TRUE_FREQS_EQUAL = [0.25, 0.25, 0.25, 0.25]
TRUE_BRANCH_LENGTHS_RAW = {"A": 0.2, "B": 0.2, "C": 0.2, "D": 0.2}

TEST_CONFIGS = [
    {
        "name": "at_mild",
        "alt_freqs_str": "0.3,0.2,0.2,0.3",
        "true_comp": [0.3, 0.2, 0.2, 0.3],
    },
    {
        "name": "at_med",
        "alt_freqs_str": "0.4,0.1,0.1,0.4",
        "true_comp": [0.4, 0.1, 0.1, 0.4],
    },
    {
        "name": "at_strong",
        "alt_freqs_str": "0.45,0.05,0.05,0.45",
        "true_comp": [0.45, 0.05, 0.05, 0.45],
    },
    {
        "name": "gc_mild",
        "alt_freqs_str": "0.2,0.3,0.3,0.2",
        "true_comp": [0.2, 0.3, 0.3, 0.2],
    },
    {
        "name": "gc_med",
        "alt_freqs_str": "0.1,0.4,0.4,0.1",
        "true_comp": [0.1, 0.4, 0.4, 0.1],
    },
    {
        "name": "gc_strong",
        "alt_freqs_str": "0.05,0.45,0.45,0.05",
        "true_comp": [0.05, 0.45, 0.45, 0.05],
    },
]

MODEL = "HKY"
TSTV = 2.0
TREE_HEIGHT = 0.5

# Raw tree: ((A:0.2,B:0.2):0.05,(C:0.2,D:0.2):0.05);
# root-to-tip = 0.25 (all four paths identical in this symmetric tree);
# with -d 0.5 → scale = 0.5 / 0.25 = 2.0
RAW_TREE_HEIGHT = 0.25
EFFECTIVE_SCALE = TREE_HEIGHT / RAW_TREE_HEIGHT  # = 2.0
LENGTHS = [200, 500, 1000, 2000, 5000]
SEEDS = [42, 43, 44]

os.makedirs(OUT_DIR, exist_ok=True)

# ==================== Helpers ====================


def run_seq_sim(length, seed, out_dir, tag, alt_freqs_str):
    """Run seq_sim_v2 to generate .nex alignment under NDCH."""
    tree_path = os.path.join(out_dir, f"{tag}_tree.nwk")
    with open(tree_path, "w") as f:
        f.write(ANNOTATED_TREE + "\n")

    out_prefix = os.path.join(out_dir, tag)
    cmd = [
        sys.executable, "-m", "seq_sim_v2",
        "-m", MODEL, "-t", str(TSTV),
        "-l", str(length), "-z", str(seed),
        "-d", str(TREE_HEIGHT),
        "--alt-freqs", alt_freqs_str,
        "-o", "n", "-y", out_prefix, "-q",
        tree_path,
    ]

    # seq_sim_v2 must be importable as a top-level package — its parent
    # directory needs to be on PYTHONPATH.
    env = os.environ.copy()
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = PROJ_PARENT + (os.pathsep + existing if existing else "")

    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(f"  seq_sim_v2 STDOUT: {result.stdout}")
        print(f"  seq_sim_v2 STDERR: {result.stderr}")
        raise RuntimeError(f"seq_sim_v2 failed with code {result.returncode}")
    nex_file = os.path.join(out_prefix, os.path.basename(out_prefix) + ".nex")
    return nex_file, tree_path


def generate_p4_script(out_dir, tag, nex_path, clean_tree_path, true_comp):
    """Generate a p4 .py script that:
      - reads data and clean tree
      - sets up identical heterogeneous model structure
      - runs optLogLike()
      - writes per-node estimated comp and blen to a JSON file
    """
    result_json = os.path.abspath(os.path.join(out_dir, f"{tag}_p4_result.json"))
    clean_tree_abs = os.path.abspath(clean_tree_path)
    nex_abs = os.path.abspath(nex_path)
    comp_init_val = str(true_comp)

    script = f'''# {tag} — p4 NDCH parameter recovery
var.warnReadNoFile = 0
var.verboseRead = 0

read(r"{nex_abs}")
read(r"{clean_tree_abs}")
t = var.trees[0]
t.taxNames = ['A','B','C','D']
t.data = Data()

# --- Determine p4 node numbering ---
root = t.root
ab_parent = None
cd_parent = None
for n in t.iterNodesNoRoot():
    leaves = []
    def _get_leaves(node, acc):
        if node.isLeaf:
            acc.append(node.name)
        else:
            c = node.leftChild
            while c:
                _get_leaves(c, acc)
                c = c.sibling
    _get_leaves(n, leaves)
    leaves_sorted = sorted(leaves)
    if leaves_sorted == ['A','B']:
        ab_parent = n
    elif leaves_sorted == ['C','D']:
        cd_parent = n

if ab_parent is None or cd_parent is None:
    import sys
    sys.exit("ERROR: could not find AB/CD clade parents in p4 tree")

print(f"AB parent node: {{ab_parent.nodeNum}}")
print(f"CD parent node: {{cd_parent.nodeNum}}")

# --- Set up heterogeneous model ---
# c0: equal freqs (fixed, f=0)
c0 = t.newComp(free=0, spec='equal')
# c1: alternative composition (free — p4 estimates this from data)
c1 = t.newComp(free=1, spec='specified', val={comp_init_val})

# Default: entire tree uses equal freqs (c0)
t.setModelComponentOnNode(c0, node=root, clade=1)
# Override: AB clade uses alternative composition (c1)
t.setModelComponentOnNode(c1, node=ab_parent, clade=1)

# Shared HKY rate matrix (fixed at kappa=2.0)
r = t.newRMatrix(free=0, spec='2p', val=2.0)
t.setNGammaCat(nGammaCat=1)
t.setPInvar(free=0, val=0.0)

# --- ML optimization ---
print("Running optLogLike ...")
t.optLogLike(verbose=1)

# --- Extract results ---
out = {{
    "logL": t.logLike,
    "comps": {{}},
    "blens": {{}},
    "free_comps": {{}},
}}

for comp in t.model.parts[0].comps:
    idx = comp.num
    out["comps"][str(idx)] = list(comp.val)

for n in t.iterNodesNoRoot():
    if n.isLeaf:
        label = n.name
    elif n.nodeNum == ab_parent.nodeNum:
        label = "AB_internal"
    elif n.nodeNum == cd_parent.nodeNum:
        label = "CD_internal"
    else:
        label = "int%d" % n.nodeNum
    out["blens"][label] = n.br.len

for n in t.iterNodes():
    if n.parts:
        out["free_comps"][str(n.nodeNum)] = {{
            "compNum": n.parts[0].compNum,
            "isleaf": n.isLeaf,
            "name": n.name if n.isLeaf else f"int_{{n.nodeNum}}"
        }}

import json
with open(r"{result_json}", "w") as f:
    json.dump(out, f, indent=2)

print(f"Results written to {result_json}")
'''
    script_path = os.path.join(out_dir, f"{tag}_p4.py")
    with open(script_path, "w") as f:
        f.write(script)
    return script_path, result_json


def run_p4(script_path):
    """Execute p4 script via the p4 binary."""
    cmd = [P4_BIN, script_path]
    env = os.environ.copy()
    p4_root = os.path.dirname(os.path.dirname(P4_BIN))  # .../bin/p4 → p4-phylogenetics-master/
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = p4_root + (os.pathsep + existing if existing else "")
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(f"P4 STDOUT: {result.stdout}")
        print(f"P4 STDERR: {result.stderr}")
        raise RuntimeError(f"p4 failed with code {result.returncode}")
    return result.stdout


def parse_p4_results(result_json):
    """Load p4 JSON output."""
    with open(result_json) as f:
        return json.load(f)


def compare_results(p4_out, tag, length, seed, true_comp):
    """Compare p4-estimated parameters with true values."""
    errors = {}

    # 1. Branch length comparison
    true_blens_eff = {k: v * EFFECTIVE_SCALE for k, v in TRUE_BRANCH_LENGTHS_RAW.items()}
    blens = p4_out.get("blens", {})
    blen_errs = {}
    for tip, true_bl in true_blens_eff.items():
        if tip in blens:
            blen_errs[tip] = abs(blens[tip] - true_bl)
    if blen_errs:
        errors["blen_mae"] = float(np.mean(list(blen_errs.values())))
        true_bls = [true_blens_eff[t] for t in blens if t in true_blens_eff]
        est_bls = [blens[t] for t in blens if t in true_blens_eff]
        if len(true_bls) >= 2 and np.std(true_bls) > 1e-10:
            errors["blen_pearson"] = float(np.corrcoef(true_bls, est_bls)[0, 1])
        else:
            errors["blen_pearson"] = float('nan')

    # 2. Composition comparison: find the free composition closest to true_comp
    comps = p4_out.get("comps", {})
    best_dist = float('inf')
    best_vals = None
    for comp_idx_str, vals in comps.items():
        vals_list = [float(x) for x in vals]
        if all(abs(v - 0.25) < 0.02 for v in vals_list):
            continue
        dist = sum((vals_list[i] - true_comp[i])**2 for i in range(4))
        if dist < best_dist:
            best_dist = dist
            best_vals = vals_list

    if best_vals:
        comp_errors = [abs(best_vals[i] - true_comp[i]) for i in range(4)]
        errors["comp_mae"] = float(np.mean(comp_errors))
        errors["comp_estimated"] = [float(x) for x in best_vals]
        errors["comp_true"] = list(true_comp)
        errors["comp_equal_ok"] = True

    errors["logL"] = p4_out.get("logL", None)
    errors["tag"] = tag
    errors["length"] = length
    errors["seed"] = seed

    return errors


# ==================== Main ====================

def main():
    all_results = []

    for config in TEST_CONFIGS:
        cfg_name = config["name"]
        alt_freqs_str = config["alt_freqs_str"]
        true_comp = config["true_comp"]

        for length in LENGTHS:
            for seed in SEEDS:
                tag = f"{cfg_name}_l{length}_s{seed}"
                tag_dir = os.path.join(OUT_DIR, tag)
                os.makedirs(tag_dir, exist_ok=True)

                print(f"\n{'='*60}")
                print(f"  [{cfg_name}] {tag}")
                print(f"{'='*60}")

                # Step 1: seq_sim_v2
                print("  [1/3] seq_sim_v2 generating data ...")
                try:
                    nex_path, _ = run_seq_sim(length, seed, tag_dir, tag, alt_freqs_str)
                except RuntimeError as e:
                    print(f"  [FAIL] seq_sim_v2 error: {e}")
                    all_results.append({
                        "config": cfg_name, "tag": tag, "length": length,
                        "seed": seed, "status": "FAIL", "error": str(e)
                    })
                    continue

                # Write clean tree (without annotations) for p4
                clean_tree_path = os.path.join(tag_dir, f"{tag}_clean.nwk")
                with open(clean_tree_path, "w") as f:
                    f.write(CLEAN_TREE + "\n")

                # Step 2: p4
                print("  [2/3] p4 inferring parameters ...")
                script_path, result_json = generate_p4_script(
                    tag_dir, tag, nex_path, clean_tree_path, true_comp
                )
                try:
                    stdout = run_p4(script_path)
                    p4_out = parse_p4_results(result_json)
                except RuntimeError as e:
                    print(f"  [FAIL] p4 error: {e}")
                    all_results.append({
                        "config": cfg_name, "tag": tag, "length": length,
                        "seed": seed, "status": "FAIL", "error": str(e)
                    })
                    continue

                # Step 3: Compare
                print("  [3/3] comparing results ...")
                errors = compare_results(p4_out, tag, length, seed, true_comp)
                errors["config"] = cfg_name
                errors["status"] = "PASS"
                all_results.append(errors)

                # Per-run summary
                print(f"    logL = {errors.get('logL', 'N/A')}")
                print(f"    blen MAE = {errors.get('blen_mae', 'N/A'):.4f}")
                print(f"    blen Pearson r = {errors.get('blen_pearson', 'N/A')}")
                print(f"    Effective blen (true={0.2*EFFECTIVE_SCALE:.1f}): " +
                      ", ".join(f"{t}={p4_out['blens'].get(t, '?'):.3f}"
                                for t in ['A','B','C','D']))
                if "comp_mae" in errors:
                    print(f"    comp MAE = {errors['comp_mae']:.4f}")
                    print(f"    comp est = {[f'{x:.3f}' for x in errors['comp_estimated']]}")
                    print(f"    comp true= {[f'{x:.3f}' for x in errors['comp_true']]}")

    # ==================== Summary ====================
    print(f"\n{'='*60}")
    print(f"  SUMMARY")
    print(f"{'='*60}")

    pass_count = sum(1 for r in all_results if r.get("status") == "PASS")
    fail_count = sum(1 for r in all_results if r.get("status") == "FAIL")

    blen_maes = [r["blen_mae"] for r in all_results if "blen_mae" in r]
    blen_pearsons = [r["blen_pearson"] for r in all_results if "blen_pearson" in r]
    comp_maes = [r["comp_mae"] for r in all_results if "comp_mae" in r]

    print(f"\n  Tests: {pass_count} PASS, {fail_count} FAIL out of {len(all_results)}")
    if blen_maes:
        print(f"  Branch length MAE:  mean={np.mean(blen_maes):.4f}, max={np.max(blen_maes):.4f}")
    if blen_pearsons:
        print(f"  Branch length r:    mean={np.mean(blen_pearsons):.4f}, min={np.min(blen_pearsons):.4f}")
    if comp_maes:
        print(f"  Composition MAE:    mean={np.mean(comp_maes):.4f}, max={np.max(comp_maes):.4f}")

    # Per-config breakdown
    for config in TEST_CONFIGS:
        cfg_name = config["name"]
        cfg_results = [r for r in all_results if r.get("config") == cfg_name]
        cfg_comp_maes = [r["comp_mae"] for r in cfg_results if "comp_mae" in r]
        if cfg_comp_maes:
            print(f"  [{cfg_name}] comp MAE mean={np.mean(cfg_comp_maes):.4f} "
                  f"(n={len(cfg_comp_maes)})")

    # Detailed table
    print(f"\n  {'Config':<8} {'Tag':<18} {'len':>5} {'seed':>4} {'logL':>12} {'blen_MAE':>9} {'blen_r':>7} {'comp_MAE':>9}")
    print(f"  {'-'*8} {'-'*18} {'-'*5} {'-'*4} {'-'*12} {'-'*9} {'-'*7} {'-'*9}")
    for r in all_results:
        cfg = r.get("config", "?")
        tag = r.get("tag", "?")
        length = r.get("length", 0)
        seed = r.get("seed", 0)
        logl = r.get("logL", "N/A")
        blen_mae = r.get("blen_mae", None)
        blen_r = r.get("blen_pearson", None)
        comp_mae = r.get("comp_mae", None)

        logl_str = f"{logl:.2f}" if isinstance(logl, (int, float)) else str(logl)
        blen_str = f"{blen_mae:.4f}" if blen_mae is not None else "N/A"
        blen_r_str = f"{blen_r:.4f}" if blen_r is not None else "N/A"
        comp_str = f"{comp_mae:.4f}" if comp_mae is not None else "N/A"

        print(f"  {cfg:<8} {tag:<18} {length:>5} {seed:>4} {logl_str:>12} {blen_str:>9} {blen_r_str:>7} {comp_str:>9}")

    print()

    # Write CSV
    csv_path = os.path.join(OUT_DIR, "results.csv")
    if all_results:
        import csv
        keys = ["config", "tag", "length", "seed", "status", "logL", "blen_mae",
                "blen_pearson", "comp_mae", "comp_estimated", "comp_true"]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(all_results)
        print(f"  CSV written to {csv_path}")

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
