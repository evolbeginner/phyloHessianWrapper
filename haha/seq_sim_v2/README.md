# seq_sim — Monte Carlo Simulator of Molecular Sequence Evolution

**Version**: 1.0.0

seq_sim evolves nucleotide or amino acid sequences along phylogenetic trees under
standard substitution models with rate heterogeneity. It supports closed-form
models (JC69, HKY, F84), eigendecomposition-based models (GTR, JTT, WAG, PAM,
BLOSUM62, LG, GENERAL), discrete/continuous Gamma rate variation, invariable sites,
codon-position rates, branch-length scaling, partitioned data, node-heterogeneous
models, and NEXUS (BEAST/MrBayes/IQ-TREE) input.

## Quick Start

```bash
# Simulate 1000 bp of DNA under HKY with ts/tv = 4.0
python -m seq_sim_v2 -m HKY -l 1000 -t 4.0 tree.nwk

# Simulate 500 amino-acid sites under WAG
python -m seq_sim_v2 -m WAG -l 500 tree.nwk

# Add gamma rate heterogeneity (4 discrete categories, shape = 0.5)
python -m seq_sim_v2 -m GTR -a 0.5 -g 4 tree.nwk

# Normalise tree height to 0.5 subs/site (requires rooted tree)
python -m seq_sim_v2 -m HKY -d 0.5 tree.nwk

# Full manual
python -m seq_sim_v2 --manual
```

## Installation

**Requirements**: Python ≥ 3.9, NumPy, SciPy

```bash
pip install numpy scipy
```

## Attribution

This Python implementation is based on the algorithms and mathematical
formulations first published with **Seq-Gen** (Rambaut & Grassly 1997),
independently rewritten using **scipy** and **numpy**.

### Contribution Boundary

| Category | Source | What seq_sim contributes |
|---|---|---|
| **HKY / F84 closed-form P(t)** | Seq-Gen (Rambaut & Grassly 1997) | Python re-expression with 2D matrix + `numpy.cumsum` |
| **GTR / AA eigendecomposition** | Mathematical standard (time-reversible Markov models) | `scipy.linalg.eigh` on symmetrised Q-matrix, single-row `set_vector` |
| **Gamma rate heterogeneity** | Yang (1994) — algorithm | Python re-implementation using `scipy.special.gammainc` + `scipy.stats.gamma.ppf` with three-method fallback cascade |
| **Simulation engine** | Seq-Gen pre-order recursive evolution | Unified mutation loop, node-heterogeneous model registry, dynamic recursion-limit guard |
| **Tree data structures** | Seq-Gen TNode/TTree design | Python `__slots__` classes, `comp_idx`/`rmat_idx` for per-node models, iterative tree-height computation |
| **Newick parser** | — | Original recursive-descent parser on string buffer with quote support, annotation `[f=N]`/`[r=N]`, nested bracket handling |
| **NEXUS preprocessing** | — | Original lightweight regex-based NEXUS reader (TRANSLATE, `[&R]`/`[&U]`, BEAST/MrBayes compatibility) |
| **Node-heterogeneous models** | Concept inspired by P4 (Foster 2004) | Original Python model-registry with `(comp_idx, rmat_idx)` keying |
| **CLI interface** | Seq-Gen conventions (`-m`, `-l`, `-a`, `-g`, etc.) | Argparse-based modular CLI, built-in comprehensive manual, rich validation |
| **Empirical AA data** | Published literature (JTT, WAG, PAM, BLOSUM, LG) | Bundled as Python dictionaries — not copyrightable |

### License

MIT License. See [LICENSE](LICENSE) file for details.

The gamma distribution functions were inspired by the Yang (1994) discrete-gamma
method as implemented in PAML. The original PAML code is distributed under the
GNU GPL v3; this project contains **no GPL-licensed code** — all probability
computations use scipy's public API.

## References

### Software & Algorithms

- **Rambaut, A. and Grassly, N. C.** (1997) Seq-Gen: An application for the
  Monte Carlo simulation of DNA sequence evolution along phylogenetic trees.
  *Bioinformatics*, 13(3): 235–238.  
  [https://github.com/rambaut/Seq-Gen](https://github.com/rambaut/Seq-Gen)

- **Yang, Z.** (1994) Maximum likelihood phylogenetic estimation from DNA
  sequences with variable rates over sites: Approximate methods.
  *Journal of Molecular Evolution*, 39: 306–314.

### Statistical Algorithms

The gamma rate discretisation employs methods first published as:

- **AS 32** — Bhattacharjee, G. P. (1970) The incomplete gamma integral.
  *Applied Statistics*, 19: 285–287.

- **AS 70** — Odeh, R. E. and Evans, J. O. (1974) The percentage points of
  the normal distribution. *Applied Statistics*, 23: 96–97.

- **AS 91** — Best, D. J. and Roberts, D. E. (1975) The percentage points
  of the χ² distribution. *Applied Statistics*, 24: 385–388.

### Node-Heterogeneous Models

The concept of assigning different model parameters to different branches
of the tree was informed by:

- **Foster, P. G.** (2004) Modeling compositional heterogeneity.
  *Systematic Biology*, 53(3): 485–495.  
  [https://github.com/pgfoster/p4-phylogenetics](https://github.com/pgfoster/p4-phylogenetics)

### Empirical Amino Acid Models

- **JTT** — Jones, D. T., Taylor, W. R., and Thornton, J. M. (1992) The rapid
  generation of mutation data matrices from protein sequences.
  *Computer Applications in the Biosciences*, 8: 275–282.

- **WAG** — Whelan, S. and Goldman, N. (2001) A general empirical model of
  protein evolution derived from multiple protein families using a
  maximum-likelihood approach. *Molecular Biology and Evolution*, 18: 691–699.

- **PAM** — Dayhoff, M. O., Schwartz, R. M., and Orcutt, B. C. (1978) A model
  of evolutionary change in proteins. *Atlas of Protein Sequence and Structure*,
  5(Suppl. 3): 345–352.

- **BLOSUM** — Henikoff, S. and Henikoff, J. G. (1992) Amino acid substitution
  matrices from protein blocks. *Proceedings of the National Academy of Sciences*,
  89: 10915–10919.

- **LG** — Le, S. Q. and Gascuel, O. (2008) An improved general amino acid
  replacement matrix. *Molecular Biology and Evolution*, 25: 1307–1320.

### Nucleotide Substitution Models

- **JC69** — Jukes, T. H. and Cantor, C. R. (1969) Evolution of protein
  molecules. *Mammalian Protein Metabolism*, 3: 21–132.

- **HKY** — Hasegawa, M., Kishino, H., and Yano, T. (1985) Dating of the
  human-ape splitting by a molecular clock of mitochondrial DNA.
  *Journal of Molecular Evolution*, 22: 160–174.

- **F84** — Felsenstein, J. (1984) Distance methods for inferring phylogenies:
  a justification. *Evolution*, 38: 16–24.

- **GTR** — Tavaré, S. (1986) Some probabilistic and statistical problems in
  the analysis of DNA sequences. *Lectures on Mathematics in the Life Sciences*,
  17: 57–86.
