#!/usr/bin/env python3

import argparse, sys

parser=argparse.ArgumentParser(description='''Add marker patterns to SNP position file
    -m marker patterns
    -s snps
    -o output''')

parser.add_argument('-m', '--markers', type=str, required=True)
parser.add_argument('-s', '--snps', type=str, required=True)
parser.add_argument('-o', '--output', type=str, required=True)

args = parser.parse_args()

markers = {}
with open(args.markers,'r') as markerfile:
    for line in markerfile:
        marker, chromosome, cm, markertype, pattern = line.rstrip().split('\t')
        markers[marker] = pattern

output = open(args.output, 'w')

header = ''
with open(args.snps,'r') as snpsfile:
    for line in snpsfile:
        if line.startswith("Species"):
            header = line.rstrip()
            header += '\tPattern'
            print(header, file=output)
            continue
        
        line = line.rstrip()
        species, cross, scaffold, position, chromosome, cm, marker, markertype, lod, loddiff = line.split('\t')
        print('{}\t{}'.format(line,markers[marker]),file=output)
