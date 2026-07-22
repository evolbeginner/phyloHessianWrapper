### v1.1.0 - 2026-07-22
- **Improved:** speed-up of `julia_bl.jl` by memorizing P(t) for subtrees during Felsenstein's pruning algo phylo lik calc
- **Fixed:** fixed a bug in `phyloHessianWrapper.rb` (since v0.6.0 `-m LG+G` will be parsed as only LG)

### v1.0.1 - 2026-07-22
- **Improved:** speed-up of `julia_bl.jl`

### v1.0.0 - 2026-07-18
- **Fixed:** improved `Readme.md`

### v0.6.0 - 2026-07-09
- **New features:** DNA subs models allowed by `--st DNA -m GTR`

### v0.5.0 - 2026-06-08
- **New features:** `--fd_scheme forward` to enable forward finite difference in gradient and STK2004 hessian calculation (1st-order derivative). Default: central.
- **New features:** UDM subs models enabled, such as `-m POISSON+UDM0004CLR --tree_add_cmd "-mdef /mnt/hd1/home/sishuo/lab-tools/phyloHessian/substitution_model/merged_nexus/UDM_clr_iqtree_merged.nex`.

### v0.4.5 - 2026-06-06
- **New features:** `paml_order_unroot.R` to automatically convert a rooted tree into unrooted (ref.tre) in the PAML way, so that there's no need to input `ref.tre` by yourself.

### v0.4.4 - 2026-05-20
- **Improved:** Hessian calculation speedup
- **Fixed:** fix a bug about thread allocation

### v0.4.3 — 2026-04-13
- **Fixed:** fix a bug of failing to load `Dir.rb`.

### v0.4.2 — 2026-04-11
- **Fixed:** bug of `+I` model fixed (in v0.4.0 and v0.4.1 it could lead to NA for lnL thus no output for the Hessian).

### v0.4.1 — 2026-03-27
- **Fixed:** bug of not uploading some dependencies in the folder `additional_scripts/`

### v0.4.0 — 2026-03-27
- **Fixed:** bug fixed for LG4M combined with `+I`.
- **Improved:** speed up for lnL calculation via optimization of the matrix exponential.
- **Improved:** speed up for Hessian calculation particularly under the 2nd-order derivative of lnL (`--hessian_type fd`).

