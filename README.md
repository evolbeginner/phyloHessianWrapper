# phyloHessianWrapper

## Installation
Run the following to check dependencies are installed and if any not installed, run the cmd to install them.
`additional_scripts/check_dependency.rb`

# How to run phyloHessian
1. Get the species tree from MCMCTree or CODEML. This can be done by running MCMCtree with `usedata = 3` in mcmctree.ctl then run `sed !4d out.BV > ref.tre`. This step is important as *phyloHessian* uses R package `ape` the correct order.

2. Run the following wrapper script. It calls IQ-Tree or PhyML to obtain the MLEs of all parameters (tree lengths and others) given the fixed tree topology, then call *phyloHessian* generate the Gradients (*g*), Hessian matrix (*H*) for the branch lengths conditional on other parameters fixed.
`ruby ~/lab-tools/dating/phyloHessian/phyloHessianWrapper.rb -s sim/alignment/combined.fas --reftree ref.tre -m LG+G --outdir ph --force --cpu 10 `

