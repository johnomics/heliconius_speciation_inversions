#!/usr/bin/env python3

import argparse
import glob
from Bio import SeqIO

parser = argparse.ArgumentParser(description='''Calculate length histogram for multiple FASTA files

    -g glob of FASTA files
    -n name
    ''')

parser.add_argument('-g', '--globfasta', type=str, required=True)
parser.add_argument('-n', '--name', type=str)
args = parser.parse_args()

seqlens = {}
for fastafile in glob.glob(args.globfasta):
    for seq_record in SeqIO.parse(fastafile, "fasta"):
        seqlen = len(seq_record)
        if seqlen not in seqlens:
            seqlens[seqlen] = 0
        seqlens[seqlen] += 1

for seqlen in sorted(seqlens):
    if args.name:
        print(args.name, "\t", end='')
    print("{}\t{}".format(seqlen, seqlens[seqlen]))