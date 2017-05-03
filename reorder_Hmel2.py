#!/usr/bin/env python3

import sys, argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class Transfer:
    def __init__(self, pool, placement, poolid, newscaffold, orientation, scaffold, start, end):
        self.pool = pool
        self.placement = placement
        self.poolid = poolid
        self.newscaffold = newscaffold
        self.orientation = orientation
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.length = self.end-self.start+1

def nstats(scaffolds, n):
    cumulative_length = 0
    genome_length = sum(scaffolds)
    n_number = 0
    for length in sorted(scaffolds, reverse=True):
        cumulative_length += length
        n_number += 1
        if cumulative_length > genome_length*(n/100):
            return length, n_number

def write_scaffold(seqrec, name, outfile):
    seqrec.id = name
    seqrec.description="length="+str(len(seqrec))
    SeqIO.write(seqrec, outfile, 'fasta')

def get_name(chrom_scaffolds, chromosome):
    chrom_scaffold_id = len(chrom_scaffolds[chromosome]) + 1
    name = "Hmel2" + format(chromosome, '02d') + format(chrom_scaffold_id, '03d') + 'o'
    return name

def write_ordered_scaffold(seqrec, chromosome, chrom_scaffolds, outfile):
    name = get_name(chrom_scaffolds, chromosome)
    write_scaffold(seqrec, name, outfile)
    chrom_scaffolds[chromosome].append(seqrec)

def write_agp_part(part, chromosome_position, name, chromosome, partlength, orientation, notlast, agpfile):
    chrom_part = (part+1) * 2 -1
    chrom_part_end = chromosome_position + partlength - 1
    print("{}\t{}\t{}\t{}\tD\t{}\t1\t{}\t{}\t".format(chromosome, chromosome_position, chrom_part_end, chrom_part, name, partlength, orientation),file=agpfile)
    if notlast:
        print("{}\t{}\t{}\t{}\tN\t100\tfragment\tno".format(chromosome, chrom_part_end + 1, chrom_part_end + 100, chrom_part + 1),file=agpfile)

def write_new_transfer(t, chrom_scaffold, chromosome, scaffold_position, chromosome_position, chrom_scaffolds, part, refinedagp, transferfile):
    write_agp_part(part, chromosome_position, t.newscaffold, chromosome, len(chrom_scaffold), t.orientation, part != len(transfers[chromosome])-1, refinedagp )
    print('\t'.join([str(x) for x in [chromosome, chromosome_position, chromosome_position + len(chrom_scaffold) - 1, t.pool, t.placement, t.poolid, \
                     t.scaffold, t.start, t.end, t.orientation, t.newscaffold, len(chrom_scaffold), \
                     get_name(chrom_scaffolds, chromosome), scaffold_position, scaffold_position + len(chrom_scaffold) - 1 \
         ]]), file=transferfile)

parser=argparse.ArgumentParser(description='''Reorder Hmel2 and calculate genome stats based on revised scaffold map
    -s scaffolds
    -f fasta
    -o output
    ''')

parser.add_argument('-s', '--scaffolds', type=str, required=True)
parser.add_argument('-f', '--fasta', type=str)
parser.add_argument('-o', '--output', type=str)

args = parser.parse_args()

if args.fasta and not args.output or args.output and not args.fasta:
    print("Please specify both -f for input FASTA filename and -o for output FASTA filename")
    sys.exit()


# Load reference genome
if args.fasta:
    genome = SeqIO.to_dict(SeqIO.parse(open(args.fasta, 'r'),"fasta"))

# Load transitions file
transfers = defaultdict(list)
with open(args.scaffolds,'r') as scf:
    for line in scf:
        chromosome, pool, placement, poolid, newscaffold, orientation, scaffold, start, end, *rest = line.rstrip().split('\t')
        if chromosome == "Chromosome":
            continue
        transfers[int(chromosome)].append(Transfer(pool, placement, poolid, newscaffold, orientation, scaffold, start, end))

gap = SeqRecord(Seq('N' * 100))

chromlist = sorted(transfers)
# Move unplaced scaffolds to the end of the list
if chromlist[0] == 0:
    del chromlist[0]
    chromlist.append(0)

# Output refined scaffolds

orderfile = open(args.output + "_ordered.fa", 'w')
refinedfile = open(args.output + "_refined.fa", 'w')
refinedagp = open(args.output + "_refined.chromosome.agp", 'w')
transferfile = open(args.output + ".transitions.tsv", 'w')
print("Chromosome\tChromStart\tChromEnd\tPool\tPoolType\tPoolID\tHmel2Scaffold\tHmel2Start\tHmel2End\tOrientation\tHmel2Refined\tLength\tHmel2Ordered\tHmel2oStart\tHmel2oEnd", file=transferfile)

all_scaffolds = []

chrom_scaffolds = defaultdict(list)
for chromosome in chromlist:
    ok_chrom_scaffold = SeqRecord(Seq(''))
    scaffold_position = 1
    chromosome_position = 1
    for part, t in enumerate(transfers[chromosome]):
        chrom_scaffold = genome[t.scaffold][(t.start-1):t.end]
        write_scaffold(chrom_scaffold, t.newscaffold, refinedfile)
        all_scaffolds.append(t.length)

        if t.orientation is '-':
            chrom_scaffold = chrom_scaffold.reverse_complement()

        if t.placement == "ok":
            if ok_chrom_scaffold.seq:
                ok_chrom_scaffold += gap
                scaffold_position += len(gap)
                chromosome_position += len(gap)
            write_new_transfer(t, chrom_scaffold, chromosome, scaffold_position, chromosome_position, chrom_scaffolds, part, refinedagp, transferfile)
            
            ok_chrom_scaffold += chrom_scaffold
            scaffold_position += len(chrom_scaffold)
            chromosome_position += len(chrom_scaffold)
        else:
            if ok_chrom_scaffold.seq: # Output previous merged scaffold
                write_ordered_scaffold(ok_chrom_scaffold, chromosome, chrom_scaffolds, orderfile)
                chromosome_position += len(gap)
                ok_chrom_scaffold = SeqRecord(Seq(''))
                scaffold_position = 1
            write_new_transfer(t, chrom_scaffold, chromosome, scaffold_position, chromosome_position, chrom_scaffolds, part, refinedagp, transferfile)
            write_ordered_scaffold(chrom_scaffold, chromosome, chrom_scaffolds, orderfile)
            chromosome_position += len(chrom_scaffold) + len(gap)
            scaffold_position = 1

    if ok_chrom_scaffold.seq: # Output final scaffold
        write_ordered_scaffold(ok_chrom_scaffold, chromosome, chrom_scaffolds, orderfile)

refinedfile.close()
orderfile.close()

new_scaffolds = []
complete = []
orderagp = open(args.output + "_ordered.chromosome.agp", 'w')
for chromosome in chromlist:
    if len(chrom_scaffolds[chromosome]) == 1:
        complete.append(chromosome)
    chromosome_position = 1
    for part, chrom_scaffold in enumerate(chrom_scaffolds[chromosome]):
        write_agp_part(part, chromosome_position, chrom_scaffold.id, chromosome, len(chrom_scaffold), '+', part != len(chrom_scaffolds[chromosome])-1, orderagp)
        
        chromosome_position += len(chrom_scaffold) + 100
        new_scaffolds.append(len(chrom_scaffold))
orderagp.close()


print("Found {} refined scaffolds, length {}".format(len(all_scaffolds), sum(all_scaffolds)))
print("Found {} ordered scaffolds, length {}".format(len(new_scaffolds), sum(new_scaffolds)))
print("{} chromosomes complete: {}".format(len(complete), ','.join([str(x) for x in sorted(complete)])))
for scaffolds in all_scaffolds, new_scaffolds:
    for n in range(5, 100, 5):
        (nlen, nnum) = nstats(scaffolds, n)
        print("{}\t{}\t{}".format(n, nlen, nnum))
