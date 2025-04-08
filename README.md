# phyloHessianWrapper

## Installation
Run the following to check dependencies are installed and if any not installed, run the cmd to install them.
`additional_scripts/check_dependency.rb`

# How to run phyloHessian
1. Get the species tree from MCMCTree or CODEML. This can be done by running MCMCtree with `usedata = 3` in mcmctree.ctl then run `sed !4d out.BV > ref.tre`. This step is important as *phyloHessian* uses R package `ape` the correct order.

2. 
`ruby ~/lab-tools/dating/phyloHessian/phyloHessianWrapper.rb -s sim/alignment/combined.fas --outdir ph --force --reftree ref.tre --cpu 10 -m LG+G`

