#!/usr/bin/env python3

import os, sys
import argparse
from collections import defaultdict

parser=argparse.ArgumentParser(description='''
    -i input''')

parser.add_argument('-i', '--input', type=str)
args = parser.parse_args()

Inputstream = sys.stdin
if args.input:
    Inputstream = open(args.input)
    
depths=defaultdict(int)
for base in Inputstream:
    f = base.rstrip().split("\t")
    depths[int(f[2])] += 1

for depth in sorted(depths):
    print("{}\t{}".format(depth, depths[depth]))
