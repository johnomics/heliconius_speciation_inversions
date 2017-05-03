#!/usr/bin/env python3

import sys, argparse
from collections import defaultdict
from intervaltree import Interval, IntervalTree

parser=argparse.ArgumentParser(description='''Make master chromosome AGP including Hmel2 scaffold names and gaps
    -t transitions
    -a agp
    ''')

parser.add_argument('-t', '--transitions', type=str, required=True)
parser.add_argument('-a', '--agp', type=str, required=True)

args = parser.parse_args()

agp = defaultdict(IntervalTree)
for line in open(args.agp):
    f = line.rstrip().split('\t')
    scaffold, start, end = f[0], int(f[1]), int(f[2])
    agp[scaffold][start:end+1] = f

# Load transitions file
transfers = defaultdict(list)
curchrom = -1
for i, line in enumerate(open(args.transitions)):
    if i==0:
        hf = line.rstrip().split('\t')
        continue
    f = dict(zip(hf, line.rstrip().split('\t')))

    chromosome, hmel2scaffold, hmel2start, hmel2end = int(f['Chromosome']), f['Hmel2Scaffold'], int(f['Hmel2Start']), int(f['Hmel2End'])
    
    if chromosome == 0:
        continue

    if curchrom != chromosome:
        curchrom = chromosome
        chrompart = 1

    if chrompart > 1:
        print('{}\t{}\t{}\t{}\tN\t100\tfragment\tyes'.format(curchrom, chrompos, chrompos+99, chrompart))
        chrompart += 1
    
    agp_overlaps = agp[hmel2scaffold].search(hmel2start, hmel2end)
    reverse_agp = f['Orientation'] == '-'
    
    chrompos = int(f['ChromStart'])
    for i in sorted(agp_overlaps, reverse=reverse_agp):
        partstart = max(hmel2start, int(i.data[1]))
        partend   = min(hmel2end,   int(i.data[2]))
        partlength = partend - partstart + 1
        if i.data[4] == 'W':
            i.data[6] = str(partstart)
            i.data[7] = str(partend)
            if reverse_agp:
                i.data[8] = '-'

        print('{}\t{}\t{}\t{}\t{}'.format(f['Chromosome'], chrompos, chrompos + partlength - 1, chrompart, '\t'.join(i.data[4:])))
        chrompos += partlength
        chrompart += 1
