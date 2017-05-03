#!/usr/bin/env python3

import argparse, sys
from intervaltree import Interval, IntervalTree
from collections import defaultdict
from random import randrange

chromosome_lengths = [0, 17206585, 9045316, 10541528, 9662098, 9908586, 14054175, 14308859, 9320449, 8708747, 17965481, 11759272, 16327298, 18127314, 9174305, 10235750, 10083215, 14773299, 16803890, 16399344, 14871695, 13359691]


parser=argparse.ArgumentParser(description='''Calculate probability of missing inversions of various sizes
    -g gaps
    -o output
    ''')

parser.add_argument('-g', '--gaps', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)
args = parser.parse_args()


gaps = defaultdict(lambda:defaultdict(IntervalTree))

for line in open(args.gaps):
    if line.startswith('Species'):
        continue

    f = line.rstrip().split('\t')
    species, chromosome, gapstart, gapend, gaplength = f[0], int(f[1]), int(f[2]), int(f[3]), int(f[4])
    gaps[species][chromosome][gapstart:(gapend+1)] = gaplength

gapsim = defaultdict(lambda: defaultdict(int))
for i in range(10000, 1500001, 10000):
    for species in sorted(gaps):
        print(species, i, file=sys.stderr)
        for j in range(0, 1000):
            randchrom = randrange(1,22)
            randstart = randrange(1, chromosome_lengths[randchrom] - i)
            randend = randstart + i
            overlap = False
            for gap in gaps[species][randchrom].search(randstart, randend):
                if randstart < gap.begin and randend > gap.begin or randstart < gap.end and randend > gap.end:
                    overlap = True
                    break
            if overlap:
                gapsim[i][species] += 1

out = open(args.output, 'w')
print("Size\tSpecies\tProbability", file=out)
for i in sorted(gapsim):
    for species in sorted(gapsim[i]):
        print('{}\t{}\t{}'.format(i, species, gapsim[i][species]/1000), file=out)
        
out.close()