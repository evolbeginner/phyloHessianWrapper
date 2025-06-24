#! /bin/bash


function cn(){
	b=`basename $1`
	if [ ! -z "$2" -a "$2" == true ]; then
		c=${b%%.*}
	else
		c=${b%.*}
	fi
	echo $c
}


function mkdir_with_force(){
	[ -z $1 ] && return 1
	outdir=$1
	is_force=$2
	if [ -d $outdir ]; then
		if [ "$is_force" == true ]; then
			rm -rf $outdir
		else
			echo "outdir $outdir has existed! Exiting ......" >&2
			return 1
		fi
	fi
	mkdir -p $outdir
}

