#!/usr/bin/env python3

import argparse, sys, random
import multiprocessing as mp
from termcolor import colored
from collections import defaultdict, namedtuple

parser=argparse.ArgumentParser(description='''Collapse within-species cross maps to single species maps
    -m map
    -l LOD threshold
    -o outputstub
    -t threads
    -p process
    -e errors
    -c chromosome
    -s sizes
    ''')

hybrids = { 'C10':'C6', 'C11':'C6', 'C14':'C6', 'C18':'C6', 'C20':'C6', 'C23':'C6', 'C24':'C6', 'C8':'C6',
            'C32':'C26', 'C36':'C26', 'C37':'C26', 'C38':'C26', 'C51':'C26',
            'C70':'C29', 'C71':'C29', 'C72':'C29', 'C94':'C29',
            'C3':'C1'}

parser.add_argument('-m', '--map', type=str, required=True)
parser.add_argument('-l', '--lodthreshold', type=int)
parser.add_argument('-o', '--outputstub', type=str, required=True)
parser.add_argument('-t', '--threads', type=int, default=4)
parser.add_argument('-p', '--process', action='store_true')
parser.add_argument('-e', '--errors', type=str)
parser.add_argument('-c', '--chromosome', type=int)

args = parser.parse_args()

Marker = namedtuple('Marker',['chromosome','cm','pattern'])

markers = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
crosses = defaultdict(int)
crossmarkers = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
crosslines = defaultdict(lambda: defaultdict(lambda:defaultdict(str)))
chrompos_scfpos = defaultdict(tuple)

print("Loading errors")
errors = {}
if args.errors:
    with open(args.errors, 'r') as errfile:
        for line in errfile:
            cross, chromosome, chromposition = line.rstrip().split('\t')
            errors[(cross, chromosome, chromposition)] = 1

print("Loading map")
header = ''

offspring = {}
with open(args.map,'r') as mapfile:
    for i, line in enumerate(mapfile):
        if i == 0:
            header = line
            hf = line.rstrip().split('\t')
            continue

        f = dict(zip(hf, line.rstrip().split('\t')))
        if args.lodthreshold and f['LOD'] != 'NA' and float(f['LOD']) < args.lodthreshold:
            continue
        if (f['Cross'], f['Chromosome'], f['ChromPosition']) in errors:
            continue
        cross, chromosome, chromposition, markertype, cm, pattern = f['Cross'], int(f['Chromosome']), int(f['ChromPosition']), int(f['MarkerType']), f['cM'], f['Pattern']
        crosses[cross] = f['Species']
        offspring[cross] = len(pattern)
        if markertype is 2 or chromosome == 0 or chromosome != int(f['LinkageGroup']):
            continue
        markers[chromosome][chromposition][cross] = Marker(chromosome, cm, pattern)
        crossmarkers[cross][chromosome][chromposition] = Marker(chromosome, cm, pattern)
        crosslines[cross][chromosome][chromposition] = line
        chrompos_scfpos[(chromosome,chromposition)] = (f['Scaffold'], f['Position'])

crosslist = sorted(crosses)

class cMBlock:
    def __init__(self, cm, pos):
        self.cm = cm
        self.pos = [pos]
    
    def __repr__(self):
        return "{}\t{}\t{}".format(self.cm, sorted(self.pos), max(self.pos)-min(self.pos)+1)

    def extend(self, pos):
        self.pos.append(pos)

any_incongruities = False
for chromosome in sorted(markers):
    for cross in crosslist:
        incongruity = True
        while incongruity:
            incongruity = False

            poslist = sorted(crossmarkers[cross][chromosome])
            cms = []
            curcm = -1
            for i in range(len(poslist)):
                pos = poslist[i]
                this = crossmarkers[cross][chromosome][pos]
                cm = float(this.cm)
                if curcm == cm:
                    cms[-1].extend(pos)
                else:
                    curcm = cm
                    cms.append(cMBlock(cm, pos))

            deleted = []
            incongruities = []
            for i in range(0, len(cms)-1):
                if cms[i+1].cm < cms[i].cm:
                    incongruity = True
                    incongruities.append("Incongruity: {}\t{}\t{}\t{}".format(cross, chromosome, cms[i].cm, cms[i+1].cm))
                    incongruities.append('\t' + repr(cms[i]))
                    incongruities.append('\t' + repr(cms[i+1]))
                    if max(cms[i].pos) - min(cms[i].pos) < 100:
                        deleted.append(i)
                    if max(cms[i+1].pos) - min(cms[i+1].pos) < 100:
                        deleted.append(i+1)
            if not deleted:
                for i in incongruities:
                    any_incongruities = True
                    print(i)
                break
            deleted = set(deleted)
            for d in sorted(deleted, reverse=True):
                for pos in cms[d].pos:
                    del crossmarkers[cross][chromosome][pos]
                    del crosslines[cross][chromosome][pos]
                    del markers[chromosome][pos][cross]

editedtsv = open(args.outputstub + '.crosses.map.tsv', 'w')
print(header, end='', file=editedtsv)
for cross in sorted(crosslines):
    for chromosome in sorted(crosslines[cross]):
        for pos in sorted(crosslines[cross][chromosome]):
            print(crosslines[cross][chromosome][pos], end='', file=editedtsv)
editedtsv.close()

if any_incongruities and not args.process:
    print("Incongruities found, will not collapse maps (use -p to override)")
    sys.exit()



def hamming(prevp, nextp):
    if len(prevp) != len(nextp):
        print("Marker lengths are not the same!\n{}\n{}".format(prevp, nextp))
        sys.exit()
    recs = 0
    for i in range(len(prevp)):
        if prevp[i] != nextp[i]:
            recs += 1
    return(recs)

o_gts = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

def get_recombination_regions(chromosome):
    recregions = defaultdict(list)
    for cross in crosslist:
        species = crosses[cross]
        if cross not in o_gts[chromosome]:
            crosschrom = crossmarkers[cross][chromosome]
            cpos = sorted(crosschrom)
            for o in range(offspring[cross]):
                o_gts[chromosome][cross][o] = [(cp, crosschrom[cp].pattern[o]) for cp in cpos if crosschrom[cp].pattern[o] != '-']

        for o in range(offspring[cross]):
            gts = o_gts[chromosome][cross][o]
            for i in range(len(gts)-1):
                if gts[i][1] != gts[i+1][1]:
                    offspring_name = cross + '_' + str(o)
                    recregions[species].append((cross,  offspring_name, gts[i][0], gts[i+1][0]))

    recout = []
    for species in sorted(recregions):
        for region in sorted(recregions[species], key=lambda x: x[0]):
            cross = region[0]
            if cross in hybrids:
                cross = hybrids[cross]
            recout.append("{}\t{}\t{}\t{}\t{}\t{}".format(species, cross, region[1], chromosome, region[2], region[3]))
    return(recout)

pool = mp.Pool(processes=args.threads)

chromlist = sorted(markers)
if args.chromosome:
    chromlist = [args.chromosome]

print("Generating actual recombinations")
rectsv = open(args.outputstub + '.recombinations.tsv','w')
print('Species\tCross\tOffspring\tChromosome\tStart\tEnd',file=rectsv)
results = [pool.apply_async(get_recombination_regions, args=(chromosome,)) for chromosome in chromlist]
output = [p.get() for p in results]
for r in output:
    for line in r:
        print(line,file=rectsv)
