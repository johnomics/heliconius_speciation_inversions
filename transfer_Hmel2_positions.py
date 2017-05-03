#!/usr/bin/env python3

import sys, argparse, gzip
from collections import defaultdict
from intervaltree import Interval, IntervalTree

def convert_position(scaffold, position, hmel2):

    
    if scaffold not in hmel2:
        return None, None

    position = int(position)
    hmel2overlap = hmel2[scaffold].search(position)
    if not hmel2overlap:
        return None, None

    hmel2overlap = list(hmel2overlap)[0]
    
    chromosome, chromstart, chromend, orientation = hmel2overlap.data

    if orientation == '+':
        chrompos = position - hmel2overlap.begin + chromstart
    else:
        chrompos = hmel2overlap.end-1 - position + chromstart
    
    return chromosome, chrompos

parser=argparse.ArgumentParser(description='''
    -i input
    -t transitions
    -n numsites
    -r refined
    -v vcf
    ''')

parser.add_argument('-i', '--input', type=str)
parser.add_argument('-t', '--transitions', type=str, required=True)
parser.add_argument('-n', '--numsites', type=int)
parser.add_argument('-r', '--refined', action='store_true', default=False)
parser.add_argument('-v', '--vcf', action='store_true', default=False)

args = parser.parse_args()

hmel2 = defaultdict(IntervalTree)

chromlengths = defaultdict(int)
for line in open(args.transitions):
    if line.startswith("Chromosome"):
        continue
    f = line.rstrip().split('\t')
    chromosome, chromstart, chromend, pool,       pooltype, poolid,    hmel2scaffold, hmel2start, hmel2end,  orientation, hmel2refined, length,     hmel2ordered, hmel2ostart, hmel2oend = \
    int(f[0]),  int(f[1]),  int(f[2]), int(f[3]), f[4],     int(f[5]), f[6],          int(f[7]),  int(f[8]), f[9],        f[10],        int(f[11]), f[12],        int(f[13]),  int(f[14])

    if chromosome == 0:
        continue

    chromlengths[chromosome] = chromend # Will get set for every transition, but the last is the chromosome length

    if args.refined:
        scaffold, scfstart, scfend, scforientation = hmel2refined, 1, length, '+'
    else:
        scaffold, scfstart, scfend, scforientation = hmel2scaffold, hmel2start, hmel2end, orientation
    hmel2[scaffold][scfstart:(scfend+1)] = (chromosome, chromstart, chromend, scforientation)


Instream = sys.stdin
if args.input:
    if args.input[-3:] == ".gz":
        Instream = gzip.open(args.input, 'rt')
    else:
        Instream = open(args.input)

sites = 0
contigs_output = False

if args.vcf:
    vcflines = defaultdict(lambda:defaultdict(str))
    for line in Instream:
        if line.startswith('##contig'):
            if not contigs_output:
                for c in sorted(chromlengths):
                    print('##contig=<ID=chr{},length={}>'.format(c, chromlengths[c]))
                contigs_output = True
            continue
        elif line.startswith('#'):
            print(line, end='')
            continue
     
        sites += 1
        if args.numsites and sites > args.numsites:
            break
        
        f = line.rstrip().split('\t')
        chromosome, chrompos = convert_position(f[0], f[1], hmel2)
        if chromosome is not None:
            vcflines[chromosome][chrompos] = "chr"+str(chromosome)+'\t'+str(chrompos)+'\t'+'\t'.join(f[2:])
        
    for chromosome in sorted(vcflines):
        for position in sorted(vcflines[chromosome]):
            print(vcflines[chromosome][position])

else:
    outlines = defaultdict(lambda:defaultdict(lambda:defaultdict(lambda:defaultdict(str))))
    ohf = ["Species", "Cross", "Chromosome", "ChromPosition", "Scaffold", "Position", "LinkageGroup", "cM", "Marker", "MarkerType", "LOD", "LODDiff", "Pattern"]
    for i, line in enumerate(Instream):
        if i == 0:
            hf = line.rstrip().split('\t')
            continue
        
        sites += 1
        if args.numsites and sites > args.numsites:
            break
        
        f = dict(zip(hf, line.rstrip().split('\t')))
        f['LinkageGroup'] = f['Chromosome']
        chromosome, chrompos = convert_position(f['Scaffold'], f['Position'], hmel2)
        if chromosome is None:
            chromosome = 0
            chrompos = -1
        
        f['Chromosome'] = chromosome
        f['ChromPosition'] = chrompos
        outlines[f['Species']][f['Cross']][chromosome][chrompos] = '\t'.join([str(f[field]) for field in ohf])
    
    print('\t'.join(ohf))
    for species in sorted(outlines):
        for cross in sorted(outlines[species]):
            for chromosome in sorted(outlines[species][cross]):
                for pos in sorted(outlines[species][cross][chromosome]):
                    print(outlines[species][cross][chromosome][pos])
