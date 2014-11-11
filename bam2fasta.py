#! /usr/bin/env python

# Usage: bam2fasta.py <inputbam> <outputdir>

import os
import sys

import pysam
from Bio import SeqIO, Seq, SeqRecord

def main(in_file, outdir):
    if not os.path.exists(outdir):
	os.makedirs(outdir)
    clusters=[]
    current_cluster=[]

    bam_file=pysam.Samfile(in_file, "rb")
    for read in bam_file:
	if len(current_cluster)==0 or current_cluster[-1].pos + current_cluster[-1].qlen > read.pos - 1000:
		current_cluster.append(read)
	else:
		clusters.append(current_cluster)
		current_cluster=[]
    if len(current_cluster)!=0:
	clusters.append(current_cluster)

    cnum=0
    for cluster in clusters:
	if len(cluster) > 50:
		outfilename="{0}/{1}.cluster{2}.fa".format(outdir, os.path.splitext(os.path.basename(in_file))[0], cnum)
		outfile=open(outfilename, 'w')
		read_num=0
		for read in cluster:
			seq = Seq.Seq(read.seq)
			if read.is_reverse:
				seq=seq.reverse_complement()
			rec=SeqRecord.SeqRecord(seq, "/read_num={0} /gb={0} /clone_end=5' /name={1}".format(read_num, read.qname), "", "")
			SeqIO.write(rec, outfile, "fasta")
			read_num+=1
		cnum+=1

if __name__ == "__main__":
    main(*sys.argv[1:])
