#!/usr/bin/env python3
"""One-click demo: generate tree + simulate 3 frequency conditions + count bases."""

import subprocess
import os
import sys

DEMO_DIR = os.path.dirname(os.path.abspath(__file__))
PYTHON_EXE = r"D:\conda\conda_envs\stimu\python.exe"
TREE_SIM_ROOT = r"C:\Users\intot\Desktop\molecular_clock\practice\seq_stimulation\reference"
SEQ_SIM_ROOT  = r"C:\Users\intot\Desktop\0811"

ENV = os.environ.copy()
ENV["PYTHONPATH"] = os.pathsep.join([TREE_SIM_ROOT, SEQ_SIM_ROOT]
    + ENV.get("PYTHONPATH", "").split(os.pathsep))

TREE_FILE = os.path.join(DEMO_DIR, "demo_tree.nwk")
CONFIG    = os.path.join(DEMO_DIR, "demo_lognormal.yaml")

def run(cmd, *, step=""):
    print(f"\n{'='*60}")
    print(f"  {step}")
    print(f"  $ {' '.join(cmd)}")
    print(f"{'='*60}")
    proc = subprocess.run(cmd, env=ENV, capture_output=False, text=True)
    if proc.returncode != 0:
        print(f"  [FAIL] exit code {proc.returncode}")
        sys.exit(proc.returncode)

def main():
    print("=" * 50)
    print("  seq_sim_v2 Live Demo")
    print("=" * 50)

    # Step 1: generate tree
    run([PYTHON_EXE, "-m", "tree_sim",
         "-n", "5", "-s", "42",
         "-c", CONFIG, "-o", TREE_FILE],
        step="Step 1/4: Generate 5-tip tree (Yule, lognormal blen)")

    print("\n  Tree content:")
    with open(TREE_FILE) as f:
        print(f"  {f.read().strip()}")

    # Step 2-4: simulate 3 frequency conditions
    # IMPORTANT: treefile must come BEFORE -f because -f uses nargs="*"
    configs = [
        ("Step 2/4: Equal frequency (A=C=G=T=0.25)",
         "out_equal", []),
        ("Step 3/4: High GC (A=0.1 C=0.4 G=0.4 T=0.1)",
         "out_gc", ["-f", "0.1", "0.4", "0.4", "0.1"]),
        ("Step 4/4: High AT (A=0.4 C=0.1 G=0.1 T=0.4)",
         "out_at", ["-f", "0.4", "0.1", "0.1", "0.4"]),
    ]

    fasta_files = []
    for step_label, subdir, extra_args in configs:
        out_dir = os.path.join(DEMO_DIR, subdir)
        os.makedirs(out_dir, exist_ok=True)
        cmd = [PYTHON_EXE, "-m", "seq_sim_v2",
               "-m", "HKY", "-l", "2000", "-d", "0.5",
               "-z", "42", "-o", "f", "-y", out_dir, "-q",
               TREE_FILE]
        if extra_args:
            cmd.extend(extra_args)
        run(cmd, step=step_label)
        fasta_file = os.path.join(out_dir, f"{subdir}_1.fasta")
        fasta_files.append(fasta_file)

    # Step 5: count base frequencies
    print(f"\n{'='*60}")
    print(f"  Base Frequency Statistics")
    print(f"{'='*60}")
    count_py = os.path.join(DEMO_DIR, "count_freqs.py")
    run([PYTHON_EXE, count_py] + fasta_files, step="counting base frequencies")

    print("\n" + "=" * 50)
    print("  Demo complete!")
    print("=" * 50)


if __name__ == "__main__":
    main()
