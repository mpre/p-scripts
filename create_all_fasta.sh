#! /usr/bin/env bash

HDIR=$1

for dirn in $(find $HDIR -type d -name 'chr*' -print | sort --numeric-sort); do
	echo ">>>>> Current chromosome $(basename $dirn)"
	for filen in $(find $dirn -type f -name '*.bam'); do
		echo "> Working on $filen"
		bam2fasta.py -b $filen -o $(basename ${HDIR/sections/fasta})/$(basename $dirn)
	done
done
