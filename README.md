# PhyloHessianWrapper

A computational pipeline for phylogenetic Hessian matrix estimation, enabling divergence time uncertainty quantification under complex substitution models.

## Table of Contents
- [Installation](#installation)
- [Quick Start](#Quick-start)
- [Detailed Usage](#detailed-usage)
- [Output Files](#output-files)

## Installation
```bash
additional_scripts/check_dependency.rb
```

### Dependencies
Ensure these tools are installed:
- Ruby (≥ 2.0)
- R (≥ 4.0) with `ape` package
- IQ-Tree2 (≥ 2.0) or PhyML (≥ 3.1)

## Quick-start
1. Preparation
Get the species tree from `MCMCTree` or `Codeml`. This can be done by running `MCMCtree` with `usedata = 3` in mcmctree.ctl then run `sed !4d out.BV > ref.tre`. This step is important as *phyloHessian* uses R package `ape` to get the order of tips which are generally different from that used by PAML. *PhyloHessian* will automatically identify any inconsistency and translate the `ape` way to the `PAML` way.

2. Run the following 
```bash
ruby ../phyloHessianWrapper.rb \
  -s sim/alignment/combined.fas \
  --reftree ref.tre \
  -m LG+G \
  --outdir LG+G_ph \
  --force \
  --cpu 10
```
`
