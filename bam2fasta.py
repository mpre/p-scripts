#! /usr/bin/env python

# Usage: bam2fasta.py <inputbam> <outputdir>

import os
import sys

import pysam
from Bio import SeqIO, Seq, SeqRecord

def main(in_file, outdir):
    if not os.path.exists(outdir):
	os.makedirs(outdir)
    out_file = "{0}/{1}.fa".format(outdir, os.path.splitext(os.path.basename(in_file))[0])
    with open(out_file, "w") as out_handle:
        # Write records from the BAM file one at a time to the output file.
        # Works lazily as BAM sequences are read so will handle large files.
        SeqIO.write(bam_to_rec(in_file), out_handle, "fasta")

def bam_to_rec(in_file):
    """Generator to convert BAM files into Biopython SeqRecords.
    """
    bam_file = pysam.Samfile(in_file, "rb")
    read_num=0
    for read in bam_file:
        seq = Seq.Seq(read.seq)
        if read.is_reverse:
            seq = seq.reverse_complement()
        rec = SeqRecord.SeqRecord(seq, "/read_num={1} /gb={0} /clone_end=5'".format(read.qname, read_num), "", "")
	read_num += 1
        yield rec

if __name__ == "__main__":
    main(*sys.argv[1:])
