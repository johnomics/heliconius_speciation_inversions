#!/usr/bin/env python3

import sys, argparse, gzip
from collections import defaultdict

erato = defaultdict(tuple)
curchrom = -1
for line in open("genome/Herato_scaffolds.tsv"):
    scaffold, length = line.rstrip().split('\t')
    length = int(length)
    if 'mt' in scaffold:
        continue
    chromosome = int(scaffold[6:8])
    if curchrom != chromosome:
        curchrom = chromosome
        chrompos = 1
    erato[scaffold] = (chromosome, chrompos, length)

    chrompos += length + 100

print("Chromosome\tMarker\tStart\tEnd\tLength")

for chrom in range(1,22):
    curmarker = -1
    outmarkers = []
    for line in open("erato/map" + str(chrom) + ".txt"):
        scaffold, position, marker = line.rstrip().split('\t')
        position = int(position.replace('*',''))
        marker = int(marker)
        
        if scaffold not in erato:
            continue

        chromosome = erato[scaffold][0]
        chrompos   = erato[scaffold][1] + position - 1
        
        if curmarker != marker:
            if curmarker != -1:
                outmarkers.append((chromosome, curmarker, start, end, end-start+1))
            curmarker = marker
            start = chrompos
        else:
            end = chrompos
    outmarkers.append((chromosome, curmarker, start, end, end-start+1))

    ignore = {}
    for i, om in enumerate(outmarkers):
        if i==0 or i==len(outmarkers)-1:
            continue
        if outmarkers[i-1][1] == outmarkers[i+1][1] and outmarkers[i][4] < outmarkers[i-1][4] and outmarkers[i][4] < outmarkers[i+1][4]:
            ignore[(outmarkers[i][1], outmarkers[i][2])] = 1
        
    for om in outmarkers:
        if (om[1], om[2]) not in ignore:
            print('{}\t{}\t{}\t{}\t{}'.format(om[0], om[1], om[2], om[3], om[4]))