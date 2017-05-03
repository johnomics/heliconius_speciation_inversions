#!/usr/bin/env python3

import sys, argparse, os
from collections import defaultdict
from random import randint

parser=argparse.ArgumentParser(description='''
    -i inversions
    -t transitions
    -s samples
    -o outputstub
    ''')

parser.add_argument('-i', '--inversions', type=str, required=True)
parser.add_argument('-t', '--transitions', type=str, required=True)
parser.add_argument('-s', '--samples', type=int, default=0)
parser.add_argument('-o', '--outputstub', type=str, default='test')

args = parser.parse_args()

chromlengths = defaultdict(int)
for line in open(args.transitions):
    if line.startswith("Chromosome"):
        continue
    f = line.rstrip().split('\t')
    chromosome, chromstart, chromend, pool,       pooltype, poolid,    hmel2scaffold, hmel2start, hmel2end,  orientation, hmel2refined, length,     hmel2ordered, hmel2ostart, hmel2oend = \
    int(f[0]),  int(f[1]),  int(f[2]), int(f[3]), f[4],     int(f[5]), f[6],          int(f[7]),  int(f[8]), f[9],        f[10],        int(f[11]), f[12],        int(f[13]),  int(f[14])

    if chromosome == 0:
        continue

    chromlengths[chromosome] = chromend

genomestarts = defaultdict(int)
genomelength = 0
for c in sorted(chromlengths):
    genomestarts[c] = genomelength + 1
    genomelength += chromlengths[c]


def checkpos(pos, chromlen):
    if pos < 0:
        pos = 1
    if pos > chromlen:
        pos = chromlen
    return pos

def get_random_position(regionlength, chromlengths, genomestarts, genomelength):
    samplechrom = 1
    randpos = genomelength
    while randpos > chromlengths[samplechrom] - regionlength: # Ensure region does not extend beyond the chromosome end
        random_genome_pos = randint(1, genomelength)
        for c in range(1, max(chromlengths)+1):
            if random_genome_pos < (genomestarts[c] + chromlengths[c]):
                samplechrom = c
                randpos = random_genome_pos - genomestarts[c] + 1
                break

    return samplechrom, randpos

inv_ids = defaultdict(int)

if not os.path.exists(args.outputstub):
    os.makedirs(args.outputstub)

for line in open(args.inversions):
    if line.startswith('GroupID'):
        continue
    f = line.rstrip().split('\t')
    hittype, chromosome, chromstart, chromend, length, group_status = f[3], int(f[6]), int(f[7]), int(f[8]), int(f[9]), f[15]

    if hittype == "Group":
        inv_ids[chromosome] += 1
        inv_name = 'chr{}.{}'.format(chromosome, inv_ids[chromosome])
        
        winout = open('{0}/{0}.{1}.tsv'.format(args.outputstub, inv_name), 'w')

        for win in range(-5,6):
            offset = win * length
            start = checkpos(chromstart + offset, chromlengths[chromosome]) # + is right here because win is negative, so offset is negative
            end = checkpos(chromend + offset, chromlengths[chromosome])
            if start == end:
                continue
            print("chr{} {} {} {} {} {}".format(chromosome, start, end, inv_name, win, 0), file=winout)

        regionlength = length*11
        for s in range(1,args.samples+1):
            samplechrom, randpos = get_random_position(regionlength, chromlengths, genomestarts, genomelength)

            for i in range(11):
                start = randpos + i*length
                print("chr{} {} {} {} {} {}".format(samplechrom, start, start+length-1, inv_name, i-5, s), file=winout)

        winout.close()
