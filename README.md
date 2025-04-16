phyloHessianWrapper.rb

# PhyloHessianWrapper

A computational pipeline for phylogenetic Hessian matrix estimation, enabling divergence time uncertainty quantification under complex substitution models.

## Table of Contents
- [Installation](#installation)
- [Quick Start](#Quick-start)
- [Detailed Usage](#detailed-usage)
- [Output Files](#output-files)
- [Literature cited](#Literature-cited)

## Installation
```bash
additional_scripts/check_dependency.rb
```

### Dependencies
Ensure these tools are installed:
- Julia (>=1.6)
- Ruby (≥ 2.7)
- R (≥ 4.0)
- IQ-Tree2 (≥ 2.0) or PhyML (≥ 3.1)

## Quick-start
1. Preparation

Get the species tree from `MCMCTree` or `Codeml`. This can be done by running `MCMCtree` with `usedata = 3` in mcmctree.ctl then run `sed !4d out.BV > ref.tre`. This step is important as *phyloHessian* uses R package `ape` to get the order of tips which are generally different from that used by PAML. *PhyloHessian* will automatically identify any inconsistency and translate the `ape` way to the `PAML` way.

In the given example, `ref.tre` can be obtained by the following bash cmd where `example/dating/ori/combined/in.BV` is the in.BV file generated by running `Codeml` or `MCMCTree` with the default settings.

```bash
sed '4!d' example/dating/ori/combined/in.BV
```

2. Run the following (with **LG+G model** and compare it with the result of `Codeml`)

Enter the folder `example/` and run *phyloHessianWrapper.rb*.

```bash
cd example/;

ruby ../phyloHessianWrapper.rb \
  -s sim/alignment/combined.fas \
  --reftree ref.tre \
  -m LG+G \
  --outdir LG+G_ph \
  --force \
  --cpu 10
```

3. Check results (optional)

Compare the lnL of the tree calculated by *phyloHessian* from `LG+G_ph/inBV/info` and that calculated by IQ-Tree (`LG+G_ph/iqtree/iqtree.log`). They should be almost the same.

4. run MCMCTree

Copy the file `LG+G_ph/inBV/in.BV` to the folder for MCMCTree `dating/LG+G/combined/`. Make sure that `usedata = 2 in.BV 3` is specified in the control file `mcmctree.ctl`. Then run in Bash
```bash
mcmctree mcmctree.ctl
```

Then compare the results with those in the folder `ori/` which uses Codeml to calculate the Hessian. Don't they look very similar?

5. run MCMCTree with a more complex model say LG+C60+G+I by specifying `-m LG+G` for #phyloHessianWrapper.rb#.


## Literature-cited
Wang S. phyloHessian: a tool to calculate phylogenetic Hessian under various substitution models and links it to MCMCTree molecular clock-dating.

