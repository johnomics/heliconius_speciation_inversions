#!/usr/bin/env python3

import argparse, os.path
from collections import defaultdict

samples = ['cydno_females','cydno_males','melpomene_females','melpomene_males']

parser=argparse.ArgumentParser(description='''Compile and output tails from all samples
    -i inputsuffix
    -t transitions
    ''')

parser.add_argument('-i', '--inputsuffix', type=str, required=True)
parser.add_argument('-t', '--transitions', type=str, required=True)
parser.add_argument('-s', '--samplename', type=str)
args=parser.parse_args()

class Tail:
    def __init__(self, sample, cluster, scf1, bp1, q1, scf2, bp2, q2, remainSeq, svtype, reads, zmws, evidence, initial_lengths, tail_lengths):
        self.sample = sample
        self.cluster = cluster
        self.scaffold1 = scf1
        self.breakpoint1 = int(bp1)
        self.quality1 = int(q1)
        self.scaffold2 = scf2
        self.breakpoint2 = int(bp2)
        self.quality2 = int(q2)
        self.remainSeq = int(remainSeq)
        self.svtype = svtype
        self.reads = int(reads)
        self.zmws = int(zmws)
        self.evidence = evidence
        self.initial_lengths = initial_lengths
        self.tail_lengths = tail_lengths
        
        self.overlaps = []
        self.length = self.breakpoint2 - self.breakpoint1 + 1 if self.scaffold1==self.scaffold2 else 'NA'
        self.svclass = None
        self.contained = False

scaffolds = defaultdict(tuple)
with open(args.transitions) as transfile:
    for line in transfile:
        if line.startswith("Chromosome"):
            continue
        f = line.rstrip().split('\t')
        chromosome, chromstart, chromend, pool,       pooltype, poolid,    hmel2scaffold, hmel2start, hmel2end,  orientation, hmel2refined, length,     hmel2ordered, hmel2ostart, hmel2oend = \
        int(f[0]),  int(f[1]),  int(f[2]), int(f[3]), f[4],     int(f[5]), f[6],          int(f[7]),  int(f[8]), f[9],        f[10],        int(f[11]), f[12],        int(f[13]),  int(f[14])

        if hmel2ordered not in scaffolds:
            scaffolds[hmel2ordered] = (int(chromstart), int(chromosome))

tails = []
initial_lengths=defaultdict(list)
tail_lengths=defaultdict(list)

if args.samplename:
    samples = [args.samplename]

for sample in samples:
    lengths_filename = sample + '.' + args.inputsuffix + '.verbose.lengths'
    if os.path.isfile(lengths_filename):
        for line in open(sample + '.' + args.inputsuffix + '.verbose.lengths'):
            fsample, cluster, scf1, bp1, q1, scf2, bp2, q2, remainSeq, svtype, readname, initiallength, taillength = line.rstrip().split('\t')
            initial_lengths[cluster].append(initiallength)
            tail_lengths[cluster].append(taillength)

    for line in open(sample + '.' + args.inputsuffix):
        if line.startswith('#'):
            continue
        cluster, pair, scf1, bp1, q1, scf2, bp2, q2, remainSeq, svtype, reads, zmws, evidence = line.rstrip().split('\t')
        if svtype == 'INV':
            tails.append(Tail(sample, cluster, scf1, bp1, q1, scf2, bp2, q2, remainSeq, svtype, reads, zmws, evidence, initial_lengths[cluster], tail_lengths[cluster]))

inversions = sorted([tail for tail in tails], key=lambda x: (x.scaffold1, x.breakpoint1))

inversion_df = open(args.inputsuffix + '.inversions.tsv', 'w')
print("Sample\tCluster\tScaffold1\tBreakpoint1\tMappingQuality1\tChromosome1\tChromPosition1\tScaffold2\tBreakpoint2\tMappingQuality2\tChromosome2\tChromPosition2\tRemainSeq\tReads\tLength\tInitialLengths\tTailLengths", file=inversion_df)
for i in inversions:
    chromosome1 = chromosome2 = chromposition1 = chromposition2 = 0
    if i.scaffold1 in scaffolds:
        chromosome1 = scaffolds[i.scaffold1][1]
        chromposition1 = scaffolds[i.scaffold1][0] + i.breakpoint1 - 1
    if i.scaffold2 in scaffolds:
        chromosome2 = scaffolds[i.scaffold2][1]
        chromposition2 = scaffolds[i.scaffold2][0] + i.breakpoint2 - 1

    print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(i.sample, i.cluster,
            i.scaffold1, i.breakpoint1, i.quality1, chromosome1, chromposition1,
            i.scaffold2, i.breakpoint2, i.quality2, chromosome2, chromposition2,
            i.remainSeq, i.reads, i.length, ','.join(i.initial_lengths), ','.join(i.tail_lengths)), file=inversion_df)