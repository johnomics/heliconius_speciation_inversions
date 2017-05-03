#!/usr/bin/env python3

import sys
import argparse
import glob
import subprocess
from collections import defaultdict, namedtuple

Alignment = namedtuple('Alignment', 'scaffold, position, mq, flag, orientation')

parser = argparse.ArgumentParser(description='''Calculate scaffold overlaps for Heliconius PacBio data

    -p folders of PacBio data
    -s scaffold order
    -a all
    ''')

parser.add_argument('-p', '--pacbiofolders', type=str, required=True)
parser.add_argument('-s', '--scaffoldorder', type=str, required=True)
parser.add_argument('-a', '--all', action='store_true')
args = parser.parse_args()

scaffold_length = {}
regions = open('regions.bed', 'w')
with open(args.scaffoldorder) as scaffolds:
    header = scaffolds.readline()
    for line in scaffolds:
        chromosome, pool, pooltype, poolid, newscaffold, orientation, scaffold, start, end, *more = line.rstrip().split('\t')
        start, end = int(start), int(end)
        if scaffold not in scaffold_length or end > scaffold_length[scaffold]:
            scaffold_length[scaffold] = end
        length = end-start+1
        if args.all or length < 200000:
            regions.write('{}\t{}\t{}\n'.format(scaffold, start, end))
        else:
            regions.write('{}\t{}\t{}\n'.format(scaffold, start, start+99999))
            regions.write('{}\t{}\t{}\n'.format(scaffold, end-99999, end))

regions.close()

folders = args.pacbiofolders.split(',')

reads = {}

for sample in folders:
    print(sample,file=sys.stderr)
    cur_chrom = ''
    chrom_reads = 0
    bamfile = sample + '/' + sample + '.sorted.bam'
    p = subprocess.Popen(['samtools','view','-Lregions.bed', bamfile], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

    for read in iter(p.stdout.readline, b''):
        read, flag, scaffold, position, mq, *args = read.decode('utf-8').rstrip().split('\t')
        position, mq, flag = int(position), int(mq), int(flag)
        orientation = '+'
        if flag & 16:
            orientation = '-'

        if read not in reads:
            reads[read] = {}
        if scaffold not  in reads[read]:
            reads[read][scaffold] = [] 
        reads[read][scaffold].append(Alignment(scaffold, position, mq, flag, orientation))
        chrom = scaffold[:7]
        if cur_chrom != chrom:
            if cur_chrom != '':
                print(file=sys.stderr)
            cur_chrom = chrom
            print(cur_chrom + "\t", end='', flush=True,file=sys.stderr)
            chrom_reads = 0

        chrom_reads += 1
        if chrom_reads % 1000 == 0:
            print('.', end='', flush=True,file=sys.stderr)
    print(file=sys.stderr)

overlaps = {}
for read in reads:
    if len(reads[read]) > 1:
        scaffold_overlap = ' '.join(sorted(reads[read]))
        if scaffold_overlap not in overlaps:
            overlaps[scaffold_overlap] = {}
        overlaps[scaffold_overlap][read] = []
        for scaffold in reads[read]:
            for alignment in reads[read][scaffold]:
                overlaps[scaffold_overlap][read].append(alignment)

for overlap in sorted(overlaps):
    positions = {}
    orientations = {}
    for read in overlaps[overlap]:
        read_orientations = {}
        for alignment in overlaps[overlap][read]:
            if alignment.scaffold not in positions:
                positions[alignment.scaffold] = []
            positions[alignment.scaffold].append(alignment.position)
            if alignment.scaffold not in read_orientations:
                read_orientations[alignment.scaffold] = {}
            if alignment.orientation not in read_orientations[alignment.scaffold]:
                read_orientations[alignment.scaffold][alignment.orientation] = 0
            read_orientations[alignment.scaffold][alignment.orientation] += 1
        orientation = ''
        for scaffold in sorted(read_orientations):
            if len(read_orientations[scaffold]) != 1:
                break
            orientation += list(read_orientations[scaffold].keys())[0]
        if len(orientation) == len(read_orientations):
            if orientation not in orientations:
                orientations[orientation] = 0
            orientations[orientation] += 1
            
    print('{}\t{}\t{}'.format(overlap, len(overlaps[overlap]), '\t'.join(['{}:{}'.format(x, str(orientations[x])) for x in orientations])))
    for scaffold in sorted(positions):
        pos = sorted(positions[scaffold])
        print('\t{}\t{}\t{}\t{}\t{}'.format(scaffold, scaffold_length[scaffold], str(min(pos)), str(max(pos)), ','.join([str(x) for x in pos])))