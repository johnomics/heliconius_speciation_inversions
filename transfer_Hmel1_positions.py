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
    ''')

parser.add_argument('-i', '--input', type=str)
parser.add_argument('-t', '--transitions', type=str, required=True)

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

    hmel2[hmel2scaffold][hmel2start:(hmel2end+1)] = (chromosome, chromstart, chromend, orientation)

ohf = ['Chromosome', 'ChromStart', 'ChromEnd', 'Hmel2_Scaffold', 'Hmel2_Start', 'Hmel2_End', 'Old_Scaffold', 'Old_Start', 'Old_End']
for i, line in enumerate(open(args.input)):
    if i == 0:
        hf = line.rstrip().split('\t')
        continue
    
    f = dict(zip(hf, line.rstrip().split('\t')))
    startchrom, startchrompos = convert_position(f['Hmel2_Scaffold'], f['Hmel2_Start'], hmel2)
    endchrom, endchrompos = convert_position(f['Hmel2_Scaffold'], f['Hmel2_End'], hmel2)
    if startchrom is None or endchrom is None:
        startchrom = endchrom = 0
        startchrompos = endchrompos = -1
    
    if startchrom != endchrom:
        continue
    
    f['Chromosome'] = startchrom
    f['ChromStart'] = startchrompos
    f['ChromEnd']   = endchrompos
    del f['Orientation']
    if f['ChromStart'] > f['ChromEnd']:
        f['ChromStart'], f['ChromEnd'] = f['ChromEnd'], f['ChromStart']
    print('\t'.join([str(f[field]) for field in ohf]))
