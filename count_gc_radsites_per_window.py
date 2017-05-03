#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC
from collections import defaultdict

parser=argparse.ArgumentParser(description='''Calculate GC content and RAD sites in windows across a genome
    -g genome FASTA
    -w window size
    -s step size
    -r restriction site''')

parser.add_argument('-g', '--genome', type=str, required=True)
parser.add_argument('-w', '--window', type=int, required=True)
parser.add_argument('-s', '--step', type=int, required=True)
parser.add_argument('-r', '--ressite', type=str, required=True)
args = parser.parse_args()

genome = {}
for scaffold in SeqIO.parse(args.genome, "fasta"):
    chromosome = int(scaffold.id[5:7])
    if chromosome not in genome:
        genome[chromosome] = scaffold
        genome[chromosome].id=chromosome
    else:
        genome[chromosome] += 'N'*100 + scaffold

del genome[0]

print("Chromosome\tPosition\tGC\tRessites")
for chromosome in genome:
    chromlen = len(genome[chromosome].seq)
    for pos in range(0, chromlen, args.step):
        winseq = genome[chromosome].seq[pos:(pos+args.window)]
        print('{}\t{}\t{:.1f}\t{}'.format(chromosome, pos+(args.window/2)+1, GC(winseq), winseq.count(args.ressite)))
