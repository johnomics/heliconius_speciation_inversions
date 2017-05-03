#!/usr/bin/env python3

import sys, argparse, gzip
from collections import defaultdict
from intervaltree import Interval, IntervalTree


parser=argparse.ArgumentParser(description='''
    -e erato_scaffolds
    -t transitions
    -m maf
    -o output
    ''')

parser.add_argument('-e', '--erato', type=str, required=True)
parser.add_argument('-t', '--transitions', type=str, required=True)
parser.add_argument('-m', '--maf', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)

args = parser.parse_args()

erato = defaultdict(tuple)
curchrom = -1
hera_chromlengths=defaultdict(int)
for line in open(args.erato):
    scaffold, length = line.rstrip().split('\t')
    length = int(length)
    if 'mt' in scaffold:
        continue
    chromosome = int(scaffold[6:8])
    if curchrom != chromosome:
        curchrom = chromosome
        chrompos = 1
    erato[scaffold] = (chromosome, chrompos, length)
    hera_chromlengths[chromosome] = chrompos + length - 1

    chrompos += length + 100
    
for c in sorted(hera_chromlengths):
    print(c, hera_chromlengths[c])

erafile = open("genome/Herato.chromosome_starts.tsv", 'w')
print("Scaffold\tChromosome\tLength\tChromStart", file=erafile)
for scaffold in sorted(erato):
    print('{}\t{}\t{}\t{}'.format(scaffold, erato[scaffold][0], erato[scaffold][2], erato[scaffold][1]), file=erafile)

erafile.close()

melpomene = defaultdict(tuple)

hmel_chromlengths = defaultdict(int)
for line in open(args.transitions):
    if line.startswith("Chromosome"):
        continue
    f = line.rstrip().split('\t')
    chromosome, chromstart, chromend, pool,       pooltype, poolid,    hmel2scaffold, hmel2start, hmel2end,  orientation, hmel2refined, length,     hmel2ordered, hmel2ostart, hmel2oend = \
    int(f[0]),  int(f[1]),  int(f[2]), int(f[3]), f[4],     int(f[5]), f[6],          int(f[7]),  int(f[8]), f[9],        f[10],        int(f[11]), f[12],        int(f[13]),  int(f[14])

    if chromosome == 0:
        continue

    hmel_chromlengths[chromosome] = chromend

    if hmel2ordered not in melpomene:
        melpomene[hmel2ordered] = (chromosome, chromstart, hmel2oend)
    else:
        mo = melpomene[hmel2ordered]
        melpomene[hmel2ordered] = (mo[0], mo[1], hmel2oend)

def convert_position(scaffold, start, genome):
    if scaffold not in genome:
        return None, None
    chromosome = genome[scaffold][0]
    chrompos   = genome[scaffold][1] + start - 1
    return chromosome, chrompos


outfile = open(args.output, 'w')
print("Hmel2_Chromosome\tHmel2_ChromStart\tHmel2_ChromEnd\tHmel2_Scaffold\tHmel2_ScfStart\tHmel2_ScfEnd\tHmel2_Length\tHmel2_Orientation\tHerato_Chromosome\tHerato_ChromStart\tHerato_ChromEnd\tHerato_Scaffold\tHerato_ScfStart\tHerato_ScfEnd\tHerato_Length\tHerato_Orientation", file=outfile)
outlines = defaultdict(lambda:defaultdict(str))
output = ''
for line in open(args.maf):
    if line.startswith('#'):
        continue
    if line.startswith('s'):
        s, scaffold, start, length, orientation, scflength, *args = line.rstrip().split()
        start, length, scflength = int(start), int(length), int(scflength)
        if orientation == '-':
            start = scflength - start + 1
        end = start + length

        if 'Herato' in scaffold:
            chromosome, chromstart = convert_position(scaffold, start, erato)
            if chromosome is not None:
                chromend = chromstart + length
                output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chromosome, chromstart, chromend, scaffold, start, end, length, orientation)

        if 'Hmel2' in scaffold:
            chromosome, chromstart = convert_position(scaffold, start, melpomene)
            if chromosome is not None and output is not '':
                chromend = chromstart + length
                linechrom, linepos = chromosome, chromstart
                output = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(chromosome, chromstart, chromend, scaffold, start, end, length, orientation) + output
                outlines[linechrom][linepos] = output
            output = ''

for chrom in sorted(outlines):
    for pos in sorted(outlines[chrom]):
        print(outlines[chrom][pos], file=outfile)

outfile.close()