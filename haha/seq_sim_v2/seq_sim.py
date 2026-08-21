"""seq_sim — Monte Carlo simulation of molecular sequence evolution along phylogenetic trees.
Usage: python -m seq_sim -m HKY -l 1000 -n 1 treefile"""

import argparse
import os
import sys
import time
import numpy as np

from .utils import set_seed, create_seed, random as _rand
from .core.models import setup_model, AA_MODELS, build_model_registry
from .core.gamma import rndgamma, discrete_gamma
from .core.evolve import evolve_sequences, _raise_recursion_limit
from .tree.tree_parser import parse_trees, allocate_sequences, dispose_tree
from .output.writers import write_sequences, write_ancestral, write_rates
from .output.nexus_reader import is_nexus, parse_nexus_trees

VERSION = "1.0.0"


def main():
    args = _parse_args()

    if args.seed is not None:
        set_seed(args.seed)
        seed = args.seed
    else:
        seed = create_seed()
        set_seed(seed)

    if not args.quiet:
        _print_banner(seed)

    tree_path = args.treefile
    if tree_path:
        with open(tree_path) as f:
            tree_text = f.read()
    else:
        tree_text = sys.stdin.read()

    if not tree_text.strip():
        sys.exit("Error: no tree data provided")

    rooted_flags = []
    if is_nexus(tree_text):
        translate_map, nexus_trees = parse_nexus_trees(tree_text)
        tree_text = '\n'.join(t for t, _ in nexus_trees)
        rooted_flags = [r for _, r in nexus_trees]
    else:
        translate_map = None

    model = setup_model(
        args.model,
        tstv=args.tstv,
        equal_exchangeability=(args.tstv is None),
        freqs=args.freqs,
        equal_freqs=(args.freqs is None),
        relative_rates=args.rates,
    )
    state_chars = model.state_chars
    is_nuc = model.num_states == 4

    model_registry = build_model_registry(
        args.model,
        freqs=args.freqs,
        tstv=args.tstv,
        equal_exchangeability=(args.tstv is None),
        equal_freqs=(args.freqs is None),
        relative_rates=args.rates,
        alt_freqs=args.alt_freqs_list,
        alt_tstv=args.alt_tstv,
        alt_rates=args.alt_rates_list,
    )

    rate_info = _setup_rates(args, args.total_sites)

    raw_trees = parse_trees(tree_text, translate_map)

    if rooted_flags:
        for (_, _, tree), rooted in zip(raw_trees, rooted_flags):
            if rooted is not None:
                tree.rooted = rooted

    datasets = _group_partitions(raw_trees, args.total_sites, args.max_partitions)

    # Raise recursion limit for deep trees
    max_nodes = 0
    for parts in datasets:
        for _, _, tree in parts:
            if tree.num_nodes > max_nodes:
                max_nodes = tree.num_nodes
    _raise_recursion_limit(max_nodes)

    fmt = args.format
    separate = args.separate or (fmt == 'f')
    out_prefix = args.output_prefix or "output"
    out_base = os.path.basename(out_prefix) or out_prefix
    suffix = args.suffix or _default_suffix(fmt)
    os.makedirs(out_prefix, exist_ok=True)

    insert_text = None
    if args.insert_file:
        with open(args.insert_file) as f:
            insert_text = f.read()

    total_start = time.time()
    multi_tree = len(datasets) > 1
    multi_dataset = args.num_datasets > 1

    if args.branch_rates:
        if args.normalize_branch_rates:
            _normalize_branch_rates(datasets)
    else:
        for parts in datasets:
            for _, _, tree in parts:
                for node in tree.node_list:
                    node.param = 1.0

    for ds_idx, parts in enumerate(datasets):
        tree_no = ds_idx + 1

        for part_len, part_rate, tree in parts:
            allocate_sequences(tree, part_len)

        if len(parts) > 1:
            wsum = sum(r * pl for pl, r, _ in parts)
            for i in range(len(parts)):
                parts[i] = (parts[i][0],
                           parts[i][1] * args.total_sites / wsum,
                           parts[i][2])

        if ds_idx == 0 and not args.quiet:
            _print_model_info(args, model, rate_info, fmt, out_prefix, seed)

        for rep in range(args.num_datasets):
            dataset_no = rep + 1

            if rate_info['type'] != 'none':
                _assign_categories(rate_info, args.total_sites, args.gamma,
                                   args.gamma_cats)

            k = 0
            for part_len, part_rate, tree in parts:
                scale = part_rate
                if args.tree_scale > 0:
                    if not tree.rooted:
                        raise ValueError(
                            "Tree must be rooted for -d (tree-scale) option. "
                            "Use -s (branch-scale) for unrooted trees."
                        )
                    scale *= args.tree_scale / tree.total_length
                elif args.branch_scale > 0:
                    scale *= args.branch_scale

                if not args.quiet and ds_idx == 0 and rep == 0:
                    eff_h = tree.total_length * scale
                    sys.stderr.write(
                        f"  Tree height: {tree.total_length:.3f}"
                        f"  → effective: {eff_h:.3f} subs/site\n")
                    if not args.tree_scale and tree.total_length > 2.0:
                        sys.stderr.write(
                            "  WARNING: tree height > 2.0 — sequences may be saturated.\n"
                            "  Consider using -d 0.5 to normalize tree height.\n")

                evolve_sequences(tree, k, part_len, scale, model_registry, rate_info,
                                 args.ancestor)
                k += part_len

            file_no = rep + 1 if multi_dataset else tree_no
            n_files = args.num_datasets if multi_dataset else len(datasets)

            if separate:
                padding = max(1, len(str(n_files)))
                fname = f"{out_base}_{file_no:0{padding}d}.{suffix}"
                path = os.path.join(out_prefix, fname)
                with open(path, 'w') as fv:
                    if fmt == 'n':
                        fv.write("#NEXUS\n[\nGenerated by seq_sim\n]\n\n")
                    _write_dataset(fv, args, trees=[t for _, _, t in parts],
                                   part_lens=[pl for pl, _, _ in parts],
                                   n_tips=parts[0][2].num_tips,
                                   state_chars=state_chars, is_nuc=is_nuc,
                                   multi_tree=multi_tree, tree_no=tree_no,
                                   multi_dataset=multi_dataset, dataset_no=dataset_no,
                                   insert_text=insert_text, fmt=fmt)
            else:
                fname = f"{out_base}.{suffix}"
                mode = 'w' if (ds_idx == 0 and rep == 0) else 'a'
                with open(os.path.join(out_prefix, fname), mode) as fv:
                    if fmt == 'n' and ds_idx == 0 and rep == 0:
                        fv.write("#NEXUS\n[\nGenerated by seq_sim\n]\n\n")
                    _write_dataset(fv, args, trees=[t for _, _, t in parts],
                                   part_lens=[pl for pl, _, _ in parts],
                                   n_tips=parts[0][2].num_tips,
                                   state_chars=state_chars, is_nuc=is_nuc,
                                   multi_tree=multi_tree, tree_no=tree_no,
                                   multi_dataset=multi_dataset, dataset_no=dataset_no,
                                   insert_text=insert_text, fmt=fmt)

        for _, _, tree in parts:
            dispose_tree(tree)

    if not args.quiet:
        elapsed = time.time() - total_start
        sys.stderr.write(f"Time taken: {elapsed:.3f} seconds\n")


def _write_dataset(fv, args, trees, part_lens, n_tips, state_chars, is_nuc,
                    multi_tree, tree_no, multi_dataset, dataset_no, insert_text, fmt):
    """Write one dataset's output to an already-open file handle."""
    if args.write_ancestral:
        if len(trees) > 1 and not args.separate:
            sys.stderr.write(
                "Warning: ancestral output with multiple partitions "
                "only writes first partition's sequences; use -os "
                "for separate files per dataset\n")
        write_ancestral(fv, fmt, trees[0], state_chars)
    else:
        write_sequences(fv, fmt, trees, part_lens, n_tips,
                         args.total_sites, state_chars,
                         tree_no if multi_tree else 0,
                         dataset_no if multi_dataset else 0,
                         is_nuc)

    if args.write_rates and not args.quiet and rate_info_global.get('site_rates') is not None:
        write_rates(sys.stderr, rate_info_global['site_rates'], args.total_sites)

    if insert_text:
        fv.write(insert_text + '\n')


def _print_manual():
    """Print the comprehensive usage manual to stdout and exit."""
    manual = r"""

================================================================================
  seq_sim v""" + VERSION + r""" — Monte Carlo Simulator of Molecular Sequence Evolution
================================================================================

seq_sim evolves nucleotide or amino acid sequences along a phylogenetic tree
under a specified substitution model. It supports both closed-form models
(JC69, HKY, F84) and eigendecomposition-based models (GTR, JTT, WAG, PAM,
BLOSUM62, LG, GENERAL) with rate heterogeneity (Gamma, codon-position rates,
invariable sites), branch-length scaling, partitioned data, heterogeneous
models per node, and NEXUS/BEAST/MrBayes input.

================================================================================
  I. QUICK START
================================================================================

  (1) Simulate 1000 bp of DNA under HKY with ts/tv=4.0 on a rooted tree:
        $ seq_sim -m HKY -l 1000 -t 4.0 tree.nwk

  (2) Simulate 500 amino-acid sites under the WAG model:
        $ seq_sim -m WAG -l 500 tree.nwk

  (3) Simulate with Gamma rate heterogeneity (4 discrete categories, shape=0.5):
        $ seq_sim -m GTR -a 0.5 -g 4 tree.nwk

  (4) Normalise tree height to 0.5 substitutions/site (requires rooted tree):
        $ seq_sim -m HKY -d 0.5 tree.nwk

  (5) Read tree from stdin, write FASTA to stdout:
        $ cat tree.nwk | seq_sim -m JC69 -l 2000

  (6) Use an ancestral root sequence instead of random generation:
        $ seq_sim -m HKY -l 1000 --ancestor "ACGTACGT..." tree.nwk

  (7) Simulate three independent datasets (replicates) per tree:
        $ seq_sim -m HKY -l 1000 -n 3 tree.nwk


================================================================================
  II. COMMAND-LINE PARAMETER REFERENCE
================================================================================

  POSITIONAL ARGUMENTS
    treefile              Newick or NEXUS tree file. Use "-" or omit for stdin.

  REQUIRED OPTIONS
    -m, --model NAME      Substitution model (REQUIRED). Supported values:
                          NUCLEOTIDE  — JC69, HKY, F84, GTR, REV (=GTR)
                          AMINO ACID  — JTT, WAG, PAM, BLOSUM, LG, GENERAL

  MODEL PARAMETERS
    -t, --tstv FLOAT      Transition/transversion ratio (HKY, F84 only).
                          Must be > 0. Omit for equal exchangeability.

    -f, --freqs FLOAT...  Equilibrium frequencies. 4 values for nucleotides
                          (A,C,G,T order), 20 for amino acids.
                          Omit for equal frequencies (default).

    -r, --rates FLOAT...  Exchangeability (relative rate) parameters.
                          6 values for GTR (AC,AG,AT,CG,CT,GT order) or
                          190 values for AA GENERAL model.
                          Omit for defaults (all 1.0 for nuc).

    --alt-freqs VEC       Alternative composition vector (comma-separated:
                          "0.3,0.2,0.3,0.2"). Repeat for multiple vectors.
                          Used with [f=N] annotations for per-node models.

    --alt-rates VEC       Alternative exchangeability vector (comma-
                          separated: "1,2,1,1,2,1"). Repeat for multiple.

    --alt-tstv FLOAT      Alternative ts/tv ratio. Repeat for multiple.

    --alt-freqs-file PATH File of alternative frequency vectors (one per line).

    --alt-rates-file PATH File of alternative exchangeability vectors.

  SIMULATION SIZE
    -l, --length INT      Sequence length in sites.  Default: 1000.
    -n, --num-datasets N  Number of independent replicates per tree.  Default: 1.

  RATE HETEROGENEITY (mutually exclusive groups: -c | -a/-g | none)
    -a, --gamma SHAPE     Gamma distribution shape parameter (alpha).
                          Typically 0.2–2.0. Enables continuous gamma
                          rate variation across sites.

    -g, --gamma-cats N    Number of discrete gamma categories (2–32).
                          Must be used WITH -a.  Common choices: 4, 8, 16.

    -c, --codon-rates R1 R2 R3
                          Codon-position rates (three values).
                          Automatically normalised so mean = 1.
                          Example: -c 1 0.5 2

    -i, --invariable PROP Proportion of invariable sites, 0 <= p < 1.
                          Sites with rate = 0. Can combine with -a/-g.

  BRANCH-LENGTH SCALING (mutually exclusive: -s | -d)
    -s, --branch-scale X  Multiply all branch lengths by X.  Default: 0 (off).

    -d, --tree-scale H    Scale tree so root-to-tip height = H subs/site.
                          Requires a ROOTED tree. Typical values:
                          0.1 (weak signal), 0.5 (moderate), 1.0 (saturated).
                          Tip: use this to control overall divergence level.

    -br, --branch-rates   Interpret [value] annotations in the Newick tree as
                          per-branch rate multipliers:  node.param = value.

    --normalize-branch-rates
                          Scale branch rate multipliers so their weighted
                          mean (by branch length) = 1.0.  Only makes sense
                          with -br.  Implicitly enables -br.

  OUTPUT
    -o, --format FMT      Output format:
                            p,phylip  — strict PHYLIP (10-char padded names)
                            r,relaxed — Relaxed PHYLIP (full names)
                            n,nexus   — NEXUS DATA block
                            f,fasta   — FASTA (default)
                          Append 's' for separate files per replicate:
                            ps, rs, ns, fs  (fs is always separate).

    -y, --output-prefix DIR
                          Output file/directory prefix.  Default: "output".

    -u, --suffix EXT      Output file suffix.  Defaults: .phy, .nex, .fasta

    -wa, --write-ancestral
                          Write internal (ancestral) node sequences.
                          Outputs in the same format as tips.

    -wr, --write-rates    Write per-site relative substitution rates.

    --ancestor SEQ        Use the given sequence string as the ancestral root
                          rather than generating a random root sequence.
    -x, --insert-file PATH
                          Text file to insert after each dataset in the output.

    -z, --seed INT        Random-number seed for reproducibility. Omit for
                          time-based seed (printed to stderr).

    -q, --quiet           Suppress all informational output (progress, seed,
                          model summary, timing).

  MISCELLANEOUS
    -p, --max-partitions N
                          Maximum number of partitions to concatenate into
                          one dataset.  Default: 1.

    --version             Print version and exit.
    --manual              Print this comprehensive manual and exit.
    -h, --help            Print short argument summary and exit.


================================================================================
  III. SUBSTITUTION MODELS
================================================================================

  NUCLEOTIDE MODELS (4 states: A, C, G, T)
  ─────────────────────────────────────────
  JC69    Jukes-Cantor 1969. Equal base frequencies, equal substitution rates.
          No extra parameters needed.

  HKY     Hasegawa-Kishino-Yano 1985. Closed-form transition probabilities.
          Uses equilibrium frequencies and a transition/transversion ratio.
          Parameters: -f (freqs), -t (ts/tv).  Omit -t for equal rates.

  F84     Felsenstein 1984. Similar to HKY but with a different parametrisation
          of the transition/transversion bias.
          Parameters: -f (freqs), -t (ts/tv).

  GTR/REV General Time-Reversible. Symmetric eigendecomposition of the 4x4
          instantaneous rate matrix. The most flexible nucleotide model.
          Parameters: -f (freqs), -r (6 exchangeability rates).

  AMINO ACID MODELS (20 states: ARNDCQEGHILKMFPSTWYV)
  ─────────────────────────────────────────────────────
  JTT     Jones-Taylor-Thornton 1992. Empirical model derived from large
          protein databases. Good general-purpose model.

  WAG     Whelan & Goldman 2001. Empirically estimated from globular proteins.

  PAM     Dayhoff PAM (Point Accepted Mutation) matrix. Historically important;
          note that some exchangeabilities are zero due to limited data.

  BLOSUM  BLOSUM62-based model. Widely used in protein alignment and phylogeny.

  LG      Le & Gascuel 2008. Modern empirical model, recommended for most
          protein analyses.

  GENERAL User-supplied empirical matrix. Requires:
          -f (20 frequencies) and -r (190 exchangeability rates).

  For any AA model, you may:
    - Omit -f to use the model's default frequencies.
    - Pass -f with 20 values to override frequencies.
    - Omit -r (preloaded data is used; only AAMODEL requires -r).


================================================================================
  IV. RATE HETEROGENEITY
================================================================================

  Rate heterogeneity allows different sites to evolve at different speeds.

  GAMMA DISTRIBUTION (+G)
    The relative substitution rate at each site is drawn from a Gamma
    distribution with shape = alpha and mean = 1.0.

    Continuous gamma (-a SHAPE):
      Each site receives an independent Gamma variate.
      More realistic but slower (one random draw per site per replicate).

    Discrete gamma (-a SHAPE -g CATEGORIES):
      The Gamma distribution is approximated by N equiprobable categories.
      The mean of the N rate values is exactly 1.0.
      Faster because P(t) matrices are precomputed per category.
      Typical values: -g 4 (standard), -g 8 (more precise).

    Shape parameter guidelines:
      alpha < 0.5  — strong rate variation (many slow, few fast sites)
      alpha = 1.0  — moderate variation
      alpha > 2.0  — weak variation (approaching equal rates)

  CODON-POSITION RATES (-c)
    Models three rate classes corresponding to codon positions 1, 2, 3.
    Example: -c 1 0.5 2 means position 1 = baseline, position 2 slower,
    position 3 twice as fast.  Values are auto-normalised to mean = 1.0.

    Mutually exclusive with Gamma (-a/-g).

  INVARIABLE SITES (+I)
    A fraction of sites never change (rate = 0). The remaining sites evolve
    at rates scaled by 1/(1-p_inv) to preserve the mean rate.
    Can be combined with Gamma models.
    Example: -a 0.5 -g 4 -i 0.2  → 20% invariable + Gamma(0.5, 4 cat).


================================================================================
  V. BRANCH-LENGTH SCALING
================================================================================

  BRANCH SCALE (-s)
    Multiply every branch length by a constant factor.
    -s 2.0  → branches twice as long (more divergence)
    -s 0.1  → branches one-tenth as long (very conserved)
    Use with unrooted OR rooted trees.

  TREE SCALE (-d)
    Normalise the whole tree so the root-to-tip distance (tree height)
    equals exactly H substitutions per site. Requires a rooted tree.
    All branches are scaled proportionally.

    Recommended values:
      -d 0.1  — weak signal, recent divergence
      -d 0.5  — moderate (good for most phylogenetics)
      -d 1.0  — strong signal
      -d 2.0  — saturated (many multiple hits)

    WARNING: if tree.total_length > 2.0 a saturation warning is printed.

  PER-BRANCH RATE MULTIPLIERS (-br)
    Newick trees with [value] annotations after branch lengths:
      (A:0.1[1.5],B:0.2[0.8])
    The [1.5] and [0.8] values become per-branch rate multipliers applied
    as:  effective_length = branch_length * scale * value.

    Use --normalize-branch-rates to globally rescale these multipliers
    so their weighted mean (by branch length) is exactly 1.0.


================================================================================
  VI. OUTPUT FORMATS
================================================================================

  PHYLIP (p)         Strict interleaved format. Taxon names truncated to 10
                      characters. Compatible with RAxML, PhyML, PAML.

  Relaxed PHYLIP (r) Like PHYLIP but full taxon names (space-separated).
                      Compatible with IQ-TREE, RAxML-NG.

  NEXUS (n)          NEXUS DATA block with Dimensions, Format, and Matrix
                      sections. Compatible with MrBayes, PAUP*, BEAST.

  FASTA (f)          Standard FASTA format, 72-char wrapped sequences.
                      Default format. Always writes separate files per
                      replicate when using multiple replicates.

  Separate files (s) Append 's' to format: -o ps, -o rs, -o ns, -o fs.
                      Creates one output file per replicate/dataset:
                        output_01.phy, output_02.phy, ...

  ANCESTRAL SEQUENCES (-wa)
    Writes sequences for ALL nodes (tips + internal) with auto-numbered
    labels. Internal nodes are numbered starting from (num_tips + 1).

  SITE RATES (-wr)
    Prints a table of per-site relative substitution rates to stderr.
    Site 1  0.85321
    Site 2  1.45210
    ...
    Useful for debugging or downstream analyses of rate heterogeneity.


================================================================================
  VII. NODE-HETEROGENEOUS MODELS (Advanced)
================================================================================

  seq_sim supports different substitution models on different branches
  via Newick annotations [f=N] and [r=N] placed after branch lengths:

    Composition index [f=N]:  selects the N-th frequency vector.
    Rate-matrix index [r=N]:  selects the N-th rate matrix / ts/tv set.

  Example — a tree where the outgroup uses different base frequencies:
    ((A:0.1[f=0],B:0.2[f=0]):0.05[f=0],C:0.3[f=1]):0.0;

  The default frequency vector (from -f) is index 0. Alternative vectors
  from --alt-freqs are indices 1, 2, ... The model registry maps each
  (comp_idx, rmat_idx) pair to a separate model instance.

  This feature generalises "mixture models" and "partitioned models" by
  allowing every branch to have its own substitution process.


================================================================================
  VIII. NEXUS INPUT (BEAST, MrBayes, IQ-TREE)
================================================================================

  seq_sim automatically detects NEXUS format and preprocesses:

    → Strips all [...] annotations (handles nesting).
    → Extracts the TRANSLATE table (numeric taxon IDs → names).
    → Extracts clean Newick trees from BEGIN TREES blocks.
    → Detects [&R] (rooted) and [&U] (unrooted) flags per tree.
    → Parses TREE and UTREE commands.

  Supported NEXUS sources:
    BEAST TreeAnnotator (MCC trees with [&height], [&posterior], etc.)
    MrBayes .con.tre (consensus trees)
    IQ-TREE .treefile, RAxML bipartitions (when wrapped in NEXUS)

  Case-insensitive keywords. Multiple TREES blocks are merged.


================================================================================
  IX. COMMON WORKFLOWS
================================================================================

  (a) Maximum-likelihood phylogeny simulation (standard):
        seq_sim -m GTR -a 0.5 -g 4 -l 1000 tree.nwk

  (b) Protein evolution with WAG+Gamma:
        seq_sim -m WAG -a 1.0 -g 8 -l 500 -o f protein_tree.nwk

  (c) NEXUS output for MrBayes analysis:
        seq_sim -m HKY -l 2000 -t 4.0 -o n tree.nwk

  (d) Generate multiple replicates for parametric bootstrap:
        seq_sim -m GTR -l 1000 -n 100 -o fs -z 42 tree.nwk

  (e) Set tree height to exactly 0.3 subs/site (rooted tree):
        seq_sim -m HKY -d 0.3 -l 500 rooted_tree.nwk

  (f) Simulate from BEAST posterior tree distribution:
        seq_sim -m GTR -l 1000 -a 1.0 -g 4 beast.trees

  (g) Use a fixed ancestral sequence (e.g., for hypothesis testing):
        seq_sim -m JC69 -l 100 --ancestor "$(cat ancestral.txt)" tree.nwk

  (h) Codon-model simulation with strong 3rd-position bias:
        seq_sim -m HKY -c 1 0.3 2.5 -l 999 tree.nwk
        (Note: length should be divisible by 3 for codon data.)

  (i) Heterogeneous model: two different GTR matrices on two clades:
        seq_sim -m GTR -r 1,2,1,1,2,1 --alt-rates "3,1,3,1,3,1" tree.nwk
        (Use with [r=0] and [r=1] annotations in the Newick string.)


================================================================================
  X. NOTES AND TROUBLESHOOTING
================================================================================

  1. ROOTED vs. UNROOTED TREES
     -d (tree-scale) requires a rooted tree. A rooted Newick tree has
     exactly two children at the root:  ((A,B),C);
     An unrooted tree has three:        ((A,B),C,D);
     Unrooted trees cannot be height-normalised. Use -s instead.

  2. SATURATION
     When tree height (total root-to-tip length) exceeds 2.0, sequences
     may be too saturated for phylogenetic signal. A warning is printed.
     Use -d to reduce the effective height.

  3. REPRODUCIBILITY
     Always set -z with a fixed seed for reproducible results:
       seq_sim -m HKY -l 1000 -z 42 tree.nwk
     Without -z, a time-based seed is used and printed to stderr.

  4. CODON RATE NORMALISATION
     Codon rates are automatically normalised such that the weighted mean
     over positions 1,2,3 equals 1.0:
       normalised_rate[i] = raw[i] * 3 / sum(raw)

  5. PARTITION MARKERS
     Trees may include [LENGTH] or [LENGTH,RATE] markers at the start.
     These define partition lengths and relative evolutionary rates for
     partitioned analyses (e.g., gene 1 = 300 bp, gene 2 = 700 bp).

  6. INVARIABLE + GAMMA
     When combining -i and -a/-g, the mean of the variable-site rates is
     scaled by 1/(1-p_inv) so the overall expected rate remains 1.0.

  7. LARGE FILES / PERFORMANCE
     For long sequences (>50kbp), large trees (>1000 tips), or many
     replicates, use -q to suppress progress output. The discrete gamma
     method (-g N) is significantly faster than continuous gamma.

  8. ANCESTRAL SEQUENCES
     The --ancestor string length must match -l. Characters must be valid
     state symbols for the chosen model (e.g., A,C,G,T for nucleotides;
     ARNDCQEGHILKMFPSTWYV for amino acids).

  9. NEXUS TREE NAMING
     Multiple TREES blocks in a single NEXUS file are supported and
     merged. TRANSLATE tables are also merged across blocks.

  10. MULTIPLE DATASETS / REPLICATES
      -n N generates N independent replicate datasets from each tree.
      Each replicate gets its own random root sequence and rate assignments.
      Use -o fs to write each replicate to a separate FASTA file.

================================================================================
  seq_sim v""" + VERSION + """
  https://github.com/anomalyco/opencode
================================================================================
"""
    print(manual)
    sys.exit(0)


def _parse_args():
    if '--manual' in sys.argv:
        _print_manual()

    p = argparse.ArgumentParser(
        prog="seq_sim",
        description="Monte Carlo simulation of molecular sequence evolution")
    p.add_argument("treefile", nargs="?", default=None,
                   help="Newick tree file (default: stdin)")
    p.add_argument("-m", "--model", required=True,
                   help="Substitution model: JC69, HKY, F84, GTR, JTT, WAG, PAM, BLOSUM, LG, GENERAL")
    p.add_argument("-l", "--length", type=int, default=1000,
                   help="Sequence length (default: 1000)")
    p.add_argument("-n", "--num-datasets", type=int, default=1,
                   help="Number of datasets per tree (default: 1)")
    p.add_argument("-p", "--max-partitions", type=int, default=1,
                   help="Maximum number of partitions (default: 1)")
    p.add_argument("-s", "--branch-scale", type=float, default=0.0,
                   help="Branch length scaling factor")
    p.add_argument("-d", "--tree-scale", type=float, default=0.0,
                   help="Normalize total tree height (root-to-tip) to N subs/site. "
                        "Recommended: 0.1 (weak), 0.5 (moderate), 1.0 (strong). "
                        "Requires rooted tree. Mutually exclusive with -s.")
    p.add_argument("-br", "--branch-rates", action="store_true",
                   help="Treat [value] after branch lengths as per-branch rate multipliers")
    p.add_argument("--normalize-branch-rates", action="store_true",
                   help="Scale branch rate multipliers to weighted mean = 1.0")
    p.add_argument("--alt-freqs", action="append", default=None,
                   help="Alternative composition vectors (comma-separated). Repeat for multiple. Nuc: A,C,G,T.")
    p.add_argument("--alt-rates", action="append", default=None,
                   help="Alternative exchangeability vectors (comma-separated). Repeat for multiple. GTR: 6 values AC,AG,AT,CG,CT,GT.")
    p.add_argument("--alt-tstv", type=float, action="append", default=None,
                   help="Alternative tstv ratios for HKY/F84. Repeat for multiple.")
    p.add_argument("--alt-freqs-file", type=str, default=None,
                   help="File with alternative composition vectors (one per line).")
    p.add_argument("--alt-rates-file", type=str, default=None,
                   help="File with alternative exchangeability vectors (one per line).")
    p.add_argument("-c", "--codon-rates", type=float, nargs=3, default=None,
                   help="Codon position rates: r1 r2 r3")
    p.add_argument("-a", "--gamma", type=float, default=None,
                   help="Gamma distribution shape (alpha)")
    p.add_argument("-g", "--gamma-cats", type=int, default=None,
                   help="Number of discrete gamma categories (2-32)")
    p.add_argument("-i", "--invariable", type=float, default=None,
                   help="Proportion of invariable sites (0-1)")
    p.add_argument("-f", "--freqs", type=float, nargs="*", default=None,
                   help="State frequencies (4 for nuc, 20 for AA, or 'e' for equal)")
    p.add_argument("-t", "--tstv", type=float, default=None,
                   help="Transition/transversion ratio (HKY/F84 only)")
    p.add_argument("-r", "--rates", type=float, nargs="*", default=None,
                   help="General rate matrix values (6 for GTR, 190 for GENERAL)")
    p.add_argument("-z", "--seed", type=int, default=None,
                   help="Random number seed")
    p.add_argument("-o", "--format", type=str, default="f",
                   help="Output format: p=PHYLIP, r=Relaxed PHYLIP, n=NEXUS, f=FASTA.  Add 's' for separate files")
    p.add_argument("-y", "--output-prefix", type=str, default=None,
                   help="Output file/directory prefix")
    p.add_argument("-u", "--suffix", type=str, default=None,
                   help="Output file suffix/extension")
    p.add_argument("-wa", "--write-ancestral", action="store_true",
                   help="Write ancestral sequences for internal nodes")
    p.add_argument("-wr", "--write-rates", action="store_true",
                   help="Write relative rate for each site to stderr")
    p.add_argument("--ancestor", type=str, default=None,
                   help="Ancestral root sequence (single string). "
                        "If provided, uses this as root instead of random generation.")
    p.add_argument("-x", "--insert-file", type=str, default=None,
                   help="Text file to insert after each dataset")
    p.add_argument("-q", "--quiet", action="store_true",
                   help="Suppress informational output")
    p.add_argument("--version", action="version",
                    version=f"seq_sim v{VERSION}")
    p.add_argument("--manual", action="store_true",
                    help="Print comprehensive usage manual and exit")
    args = p.parse_args()

    fmt_raw = args.format.lower().replace(' ', '')
    args.separate = False
    if fmt_raw.endswith('s'):
        args.separate = True
        fmt_raw = fmt_raw[:-1]
    if fmt_raw in ('p', 'phylip'):
        args.format = 'p'
    elif fmt_raw in ('r', 'relaxed'):
        args.format = 'r'
    elif fmt_raw in ('n', 'nexus'):
        args.format = 'n'
    elif fmt_raw in ('f', 'fasta'):
        args.format = 'f'
        args.separate = True
    else:
        p.error(f"Unknown output format: {fmt_raw}")

    args.total_sites = args.length

    if args.length <= 0:
        p.error(f"--length must be > 0, got {args.length}")

    if args.gamma_cats is not None:
        if args.gamma_cats < 2 or args.gamma_cats > 32:
            p.error("Gamma categories must be between 2 and 32")
        if args.gamma is None:
            p.error("Gamma shape (-a) required with discrete categories (-g)")

    if args.invariable is not None and (args.invariable < 0.0 or args.invariable >= 1.0):
        p.error("Invariable proportion must be >= 0 and < 1")

    if args.codon_rates and args.gamma is not None:
        p.error("Codon rates (-c) and gamma rates (-a/-g) cannot be used together")

    if args.branch_scale > 0 and args.tree_scale > 0:
        p.error("-s (branch-scale) and -d (tree-scale) cannot be used together")

    if args.codon_rates is not None:
        s = sum(args.codon_rates)
        if s <= 0:
            p.error("Codon rates must sum to a positive value")
        if abs(s - 3.0) > 1e-10:
            args.codon_rates = [r * 3.0 / s for r in args.codon_rates]

    if args.rates is not None and len(args.rates) not in (6, 190):
        p.error(f"-r/--rates expects 6 (GTR) or 190 (GENERAL AA) values, "
                f"got {len(args.rates)}")

    if args.freqs is not None and len(args.freqs) not in (4, 20):
        p.error(f"-f/--freqs expects 4 (nucleotide) or 20 (AA) values, "
                f"got {len(args.freqs)}")

    if args.normalize_branch_rates and not args.branch_rates:
        args.branch_rates = True

    args.alt_freqs_list = _parse_alt_vectors(args.alt_freqs, args.alt_freqs_file)
    args.alt_rates_list = _parse_alt_vectors(args.alt_rates, args.alt_rates_file)

    if args.alt_tstv and args.model.upper() not in ('HKY', 'F84'):
        p.error("--alt-tstv requires HKY or F84 model")
    if args.alt_rates_list and args.model.upper() not in ('GTR', 'REV') and args.model.upper() not in AA_MODELS:
        p.error(f"--alt-rates requires GTR or AA model, got {args.model}")

    return args


# Module-level storage for rate_info so _write_dataset can access it
rate_info_global = {}


def _setup_rates(args, num_sites):
    """Build the rate_info dictionary describing rate heterogeneity."""
    global rate_info_global
    info = {'type': 'none'}
    info['invariable'] = args.invariable is not None
    info['invariable_prop'] = args.invariable or 0.0

    if args.invariable is not None:
        info['invariable_sites'] = np.zeros(num_sites, dtype=np.int32)

    if args.codon_rates:
        info['type'] = 'codon'
        info['cat_rates'] = list(args.codon_rates)
        info['num_cats'] = 3
    elif args.gamma_cats:
        info['type'] = 'discrete_gamma'
        info['num_cats'] = args.gamma_cats
        info['gamma_shape'] = args.gamma
        _, info['cat_rates'] = discrete_gamma(args.gamma, args.gamma_cats)
        info['categories'] = np.zeros(num_sites, dtype=np.int32)
    elif args.gamma is not None:
        info['type'] = 'continuous_gamma'
        info['gamma_shape'] = args.gamma
        info['gamma_rates'] = np.zeros(num_sites)

    if args.write_rates:
        info['site_rates'] = np.zeros(num_sites)
    rate_info_global = info
    return info


def _assign_categories(rate_info, num_sites, gamma_shape, num_cats):
    """Assign random rate categories/gamma values for the current replicate."""
    if rate_info.get('invariable_sites') is not None:
        prop = rate_info['invariable_prop']
        inv = rate_info['invariable_sites']
        for i in range(num_sites):
            inv[i] = 1 if _rand() < prop else 0

    htype = rate_info['type']
    if htype == 'continuous_gamma':
        gr = rate_info['gamma_rates']
        for i in range(num_sites):
            gr[i] = rndgamma(gamma_shape) / gamma_shape
    elif htype == 'discrete_gamma':
        cats = rate_info['categories']
        for i in range(num_sites):
            cats[i] = int(_rand() * num_cats)


def _group_partitions(raw_trees, total_sites, max_parts):
    """Group parsed trees into datasets based on partition lengths."""
    if not raw_trees:
        return []

    datasets = []
    buf = []
    acc = 0

    for p_len, p_rate, tree in raw_trees:
        if p_len <= 0:
            buf.append((total_sites, p_rate, tree))
            datasets.append(buf)
            buf = []
            acc = 0
        else:
            buf.append((p_len, p_rate, tree))
            acc += p_len
            if acc >= total_sites:
                if acc != total_sites:
                    sys.stderr.write(
                        f"Warning: partition sum {acc} != {total_sites}\n")
                datasets.append(buf)
                buf = []
                acc = 0

    if buf:
        datasets.append(buf)
    return datasets


def _parse_alt_vectors(cli_args, file_path):
    result = []
    if cli_args:
        for s in cli_args:
            result.append([float(x) for x in s.split(',')])
    if file_path:
        with open(file_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                result.append([float(x) for x in line.replace(',', ' ').split()])
    return result or None


def _normalize_branch_rates(datasets):
    """Scale per-branch rate multipliers so weighted mean (by branch length) = 1.0."""
    total_weight = 0.0
    total_weighted_rate = 0.0
    for parts in datasets:
        for _, _, tree in parts:
            for node in tree.node_list:
                if node.branch0 is not None:
                    total_weight += node.length0
                    total_weighted_rate += node.length0 * node.param
    if total_weighted_rate <= 0:
        return
    norm = total_weight / total_weighted_rate
    for parts in datasets:
        for _, _, tree in parts:
            for node in tree.node_list:
                node.param *= norm


def _default_suffix(fmt):
    return {'p': 'phy', 'r': 'phy', 'n': 'nex', 'f': 'fasta'}.get(fmt, 'out')


def _print_banner(seed):
    sys.stderr.write(f"seq_sim v{VERSION}\n")
    sys.stderr.write(f"  Random seed: {seed}\n\n")


def _print_model_info(args, model, rate_info, fmt, prefix, seed):
    sys.stderr.write(f"Model: {args.model}\n")
    sys.stderr.write(f"  States: {model.num_states}\n")
    sys.stderr.write(f"  Sequence length: {args.total_sites}\n")
    sys.stderr.write(f"  Datasets per tree: {args.num_datasets}\n")
    sys.stderr.write(f"  Output: {prefix}/ [{fmt}]\n")

    if rate_info['type'] == 'codon':
        sys.stderr.write(f"  Codon rates: {rate_info['cat_rates']}\n")
    elif rate_info['type'] == 'discrete_gamma':
        sys.stderr.write(f"  Gamma: shape={args.gamma}, "
                         f"{args.gamma_cats} categories\n")
    elif rate_info['type'] == 'continuous_gamma':
        sys.stderr.write(f"  Gamma: shape={args.gamma} (continuous)\n")
    if rate_info['invariable']:
        sys.stderr.write(f"  Invariable sites: {args.invariable}\n")
    if args.branch_scale:
        sys.stderr.write(f"  Branch scale: {args.branch_scale}\n")
    if args.tree_scale:
        sys.stderr.write(f"  Tree scale: {args.tree_scale}\n")
    if args.branch_rates:
        sys.stderr.write(f"  Branch rates: enabled (per-branch multipliers from tree [value])\n")
    sys.stderr.write('\n')


if __name__ == "__main__":
    main()
