#! /bin/bash


##################################################
infile=$1

ref_tree=$2
if [ -z $ref_tree ]; then
  ref_tree=ref.tre
fi


##################################################
nw_topology -Ib $1 | ruby ~/lab-tools/dating/hessian/reorder_node.rb -i - --ref $ref_tree --output_branch > branch_out

Rscript ~/lab-tools/dating/hessian/align_bl.R 1.treefile branch_out > branch_out.matrix
