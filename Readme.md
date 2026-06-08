```markdown
# PhyloHessianWrapper

Estimate phylogenetic gradients/Hessian matrices for divergence-time inference with **MCMCtree** approximate likelihood ([dos Reis & Yang, 2011](#ref-dosreis_yang_2011)), including complex amino-acid models.

## Installation

```bash
ruby additional_scripts/check_dependency.rb
```

### Dependencies
- Julia (>= 1.6)
- Ruby (>= 2.7)
- R (>= 4.0)
- IQ-TREE2 (>= 2.0) or PhyML (>= 3.1)

---

## Quick Start

### 1) Prepare alignment + rooted species tree

You need:
- an alignment (FASTA)
- a rooted species tree (from your pipeline)

### 2) Convert rooted tree to unrooted reference tree (recommended)

```bash
Rscript additional_scripts/paml_order_unroot.R rooted.tre ref.tre
```

This is the preferred way to generate `ref.tre`.

**Alternative (old workflow):**
1. Run `mcmctree` first with `usedata=2` (default-style setup) to generate `in.BV`.
2. Extract line 4:
   ```bash
   sed '4!d' in.BV > ref.tre
   ```

### 3) Run `phyloHessianWrapper`

```bash
ruby phyloHessianWrapper.rb \
  -s sim/alignment/combined.fas \
  --reftree ref.tre \
  -m LG+G \
  --outdir LG+G_ph \
  --force \
  --cpu 10
```

### 4) (Optional) Validate likelihoods
Compare:
- `LG+G_ph/inBV/info` (phyloHessian)
- `LG+G_ph/iqtree/iqtree.log` (IQ-TREE)

### 5) Run MCMCtree with generated `in.BV`
Copy:
```bash
cp LG+G_ph/inBV/in.BV /path/to/your/mcmctree/run/
```

Then configure `mcmctree.ctl` to read `in.BV` (commonly `usedata = 2 in.BV 3`, depending on your PAML version), and run:
```bash
mcmctree mcmctree.ctl
```

---

## Main Arguments

### Required
- `-s <file>`: alignment in FASTA
- `--reftree <file>`: unrooted reference tree (recommended from `paml_order_unroot.R`)
- `-m <model>`: substitution model (e.g., `LG+G`, `LG+C60+G`, `EX2+G`)
- `--outdir <dir>`: output directory

### Optional
- `--hessian_type <STK2004|fd>`: Hessian estimator (default: `STK2004`)
- `--cpu <int>`: threads (default: 1)
- `--phylo_prog <iqtree|phyml>`: phylogeny engine (default: `iqtree`)
- `--pmsf`: enable PMSF approximation for mixture models
- `--no_mwopt`: disable IQ-TREE mixture-weight optimization
- `--force`: remove existing output directory

---

## Outputs

- `julia/`: tree-structure intermediate files
- `iqtree/` or `phyml/`: phylogenetic program outputs
- `bl/`: branch-length related files
- `ph/inBV/`:
  - `info`: log-likelihood at MLEs
  - `in.BV`: branch lengths, gradients, Hessian (for MCMCtree)
- `execution_{TIME}.log`: runtime log

---

## More Examples
### LG+C60+G with PMSF (recommended for large datasets)

You can enable PMSF in either of two equivalent ways:

```bash
# Option 1: include +PMSF in the model string
ruby phyloHessianWrapper.rb -s sim/alignment/combined.fas --reftree ref.tre -m LG+C60+G+PMSF --outdir outdir --force --cpu 10
```

```bash
# Option 2: use --pmsf flag
ruby phyloHessianWrapper.rb -s sim/alignment/combined.fas --reftree ref.tre -m LG+C60+G --pmsf --outdir outdir --force --cpu 10
```

### LG+C60+G without mixture-weight optimization

```bash
ruby phyloHessianWrapper.rb -s sim/alignment/combined.fas --reftree ref.tre -m LG+C60+G --outdir outdir --force --cpu 10 --no_mwopt
```

### EX2+G with PhyML

```bash
ruby phyloHessianWrapper.rb -s sim/alignment/combined.fas --reftree ref.tre -m EX2+G --outdir outdir --force --cpu 10 --phylo_prog phyml
```

---

## Literature Cited

### Core method / usage context
- Wang S, Meade A. 2026. *Molecular Clock Dating Using Complex Mixture Models: Applied to Ancient Symbionts*. **Mol Biol Evol**. msag039.

### MCMCtree approximate-likelihood reference
- <a id="ref-dosreis_yang_2011"></a>dos Reis M, Yang Z. 2011. Approximate likelihood calculation on a phylogeny for Bayesian estimation of divergence times. **Mol Biol Evol** 28(7):2161–2172.  
  https://academic.oup.com/mbe/article/28/7/2161/1051613

### Tools
- Minh BQ, Schmidt HA, Chernomor O, Schrempf D, Woodhams MD, von Haeseler A, Lanfear R. 2020. IQ-TREE 2. **Mol Biol Evol** 37(5):1530–1534.
- Revell LJ. 2024. phytools 2.0. **PeerJ** 12:e16505.

### Additional
- <a id="ref-Wang_2018"></a>Wang HC, Minh BQ, Susko E, Roger AJ. 2018. PMSF site-heterogeneity approximation. **Syst Biol** 67(2):216–235.
```
