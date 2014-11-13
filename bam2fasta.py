#! /usr/bin/env python

import os
import sys
import argparse
import pysam

from Bio import SeqIO, Seq, SeqRecord

def main( ):
    parser = argparse.ArgumentParser( prog = "bam2fasta.py",
                                      description = "Convert and cluster reads from bam",
                                      formatter_class = argparse.ArgumentDefaultsHelpFormatter )
    parser.add_argument( '-b', '--bam-file',
                         help = "BAM file",
                         required = False, dest = 'in_file',
                         default = "-")
    parser.add_argument( '-o', '--out-dir',
                         help = "Reads length",
                         required = False, dest = 'outdir' )
    parser.add_argument( '-c', '--cluster-min-size',
                         help = "Minimum cluster size",
                         required = False, dest = 'minsize', type = int,
                         default = 50 )
    parser.add_argument( '-m', '--max-linkage',
                         help = "Max distance between two different reads in the same cluster",
                         required = False, dest = 'linkage', type = int,
                         default = 2000 )

    args = parser.parse_args()
    outdir = args.outdir if args.outdir is not "" else "output"
    in_file = args.in_file
    linkage = args.linkage
    cluster_min_size = args.minsize

    if not os.path.exists( outdir ):
	os.makedirs( outdir )
    clusters = []
    current_cluster = []

    bam_file = pysam.Samfile( in_file, "rb" )
    for read in bam_file:
	if len( current_cluster ) == 0 or current_cluster[-1].pos + current_cluster[-1].qlen > read.pos - linkage:
            current_cluster.append( read )
	else:
            clusters.append( current_cluster )
            current_cluster = [read]
    if len( current_cluster ) != 0:
	clusters.append( current_cluster )

    cnum = 0
    for cluster in clusters:
	if len(cluster) > cluster_min_size:
            outfile = open( "{0}/{1}.cluster{2}.fa".format( outdir,
                                                            os.path.splitext( os.path.basename( in_file ) )[0] if in_file != "-" else "output",
                                                            cnum ),
                            'w' )
            read_num = 0
            for read in cluster:
                seq = Seq.Seq( read.seq )
                if read.is_reverse:
                    seq = seq.reverse_complement()
                rec = SeqRecord.SeqRecord( seq,
                                           "/read_num={0} /gb={0} /clone_end={2} /name={1}".format( read_num, read.qname,
                                                                                                    "5'" if read.is_reverse else "3'" ),
                                           "", "" )
                SeqIO.write( rec, outfile, "fasta" )
                read_num += 1
            cnum += 1

if __name__ == "__main__":
    main( )
