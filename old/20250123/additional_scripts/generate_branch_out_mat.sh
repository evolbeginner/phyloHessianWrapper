#! /bin/bash


##################################################
source ~/tools/self_bao_cun/packages/bash/util.sh


##################################################
ALIGN_BL=~/lab-tools/dating/hessian/align_bl.R


##################################################
treefile=''
ref_tree='ref.tre'
outdir=.
is_force=false


##################################################
while [ $# -gt 0 ]; do
	case $1 in
		-t|--treefile|-i|--infile)
			treefile=$2
			shift
			;;
		-r|--ref_tree|--reftree)
			ref_tree=$2
			shift
			;;
		--outdir)
			outdir=$2
			shift
			;;
		--force)
			is_force=true
			;;
	esac
	shift
done


##################################################
if [ $outdir != . ]; then
	mkdir_with_force $outdir $is_force
else
	:
fi


##################################################
nw_topology -Ib $treefile | ruby ~/lab-tools/dating/hessian/reorder_node.rb -i - --ref $ref_tree --output_branch > $outdir/branch_out

Rscript $ALIGN_BL $treefile $outdir/branch_out > $outdir/branch_out.matrix


