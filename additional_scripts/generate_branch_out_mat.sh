#! /bin/bash


##################################################
dir=`dirname $0` #additional_scripts/


##################################################
source $dir/util.sh


##################################################
ALIGN_BL=$dir/align_bl.R
REORDER_NODE=$dir/reorder_node.rb


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
nw_topology -Ib $treefile | ruby $REORDER_NODE -i - --ref $ref_tree --output_branch > $outdir/branch_out

echo "Rscript $ALIGN_BL $treefile $outdir/branch_out > $outdir/branch_out.matrix"
Rscript $ALIGN_BL $treefile $outdir/branch_out > $outdir/branch_out.matrix


