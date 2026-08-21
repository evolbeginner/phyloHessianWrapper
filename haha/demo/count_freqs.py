#!/usr/bin/env python3
"""Count A/C/G/T frequencies from one or more FASTA files."""

import sys
from collections import Counter

BASES = "ACGT"


def count_fasta(fpath):
    seqs = {}
    current_name = None
    current_seq = []
    with open(fpath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    seqs[current_name] = "".join(current_seq)
                current_name = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        if current_name is not None:
            seqs[current_name] = "".join(current_seq)
    return seqs


def report(fpath, seqs):
    print(f"\n{'='*60}")
    print(f"  File: {fpath}")
    print(f"{'='*60}")
    print(f"  Sequences: {len(seqs)}")
    print()

    total_counts = Counter()
    total_len = 0

    for name, seq in seqs.items():
        counts = Counter(seq.upper())
        seq_len = len(seq)
        total_len += seq_len
        for b in BASES:
            total_counts[b] += counts.get(b, 0)

        base_info = ", ".join(
            f"{b}={counts.get(b, 0):>5d} ({counts.get(b, 0)/seq_len*100:5.1f}%)"
            if seq_len > 0 else f"{b}=0"
            for b in BASES
        )
        print(f"  {name:<20s}  {base_info}")

    print(f"\n  {'-'*55}")
    print(f"  {'TOTAL':<20s}  ", end="")
    total_info = ", ".join(
        f"{b}={total_counts[b]:>5d} ({total_counts[b]/total_len*100:5.1f}%)"
        if total_len > 0 else f"{b}=0"
        for b in BASES
    )
    print(total_info)
    print(f"  Total bases: {total_len}")
    print()


def main():
    if len(sys.argv) < 2:
        print("Usage: python count_freqs.py <fasta_file> [fasta_file ...]")
        sys.exit(1)

    for fpath in sys.argv[1:]:
        seqs = count_fasta(fpath)
        report(fpath, seqs)


if __name__ == "__main__":
    main()
