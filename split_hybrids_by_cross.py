#!/usr/bin/env python3

import os
from collections import defaultdict

hybrids = { 'C10':'C6', 'C11':'C6', 'C14':'C6', 'C18':'C6', 'C20':'C6', 'C23':'C6', 'C24':'C6', 'C8':'C6',
            'C32':'C26', 'C36':'C26', 'C37':'C26', 'C38':'C26', 'C51':'C26',
            'C70':'C29', 'C71':'C29', 'C72':'C29', 'C94':'C29',
            'C3':'C1'}

with open('header.txt') as h:
    bc_header = h.readline().rstrip().split('\t')[3:]

bc_file = open('linkage_map.hybrids.backcross.hmel2.tsv', 'w')
f1_file = open('linkage_map.hybrids.f1s.hmel2.tsv', 'w')
all_file = open('linkage_map.hybrids.all.hmel2.tsv', 'w')

header = "Species\tCross\tScaffold\tPosition\tChromosome\tcM\tMarker\tMarkerType\tLOD\tLODDiff\tPattern"
for outfile in bc_file, f1_file, all_file:
    print(header, file=outfile)

def output_pattern(crosstype, this_pattern, last_pattern, crossfile, cm=None):
    for f in sorted(this_pattern):
        if crosstype is "Backcross" and '-' in this_pattern[f]:
            continue
        if f not in last_pattern:
            last_pattern[f] = this_pattern[f]
            if crosstype is "Backcross":
                cm[f] = 0
        else:
            if last_pattern[f] != this_pattern[f]:
                if crosstype is "Backcross":
                    cm[f] += sum(c1 != c2 for c1, c2 in zip(last_pattern[f], this_pattern[f])) / len(last_pattern[f])*100
                last_pattern[f] = this_pattern[f]
        output_cm = 'NA'
        if cm is not None:
            if crosstype is "Backcross":
                output_cm = '{:.3f}'.format(cm[f])
            else: # crosstype is "All"
                output_cm = '{:.3f}'.format(float(cm))
        print("Hybrid\t{}\t{}\t{}\t{}\t{}\tNA\t1\tNA\tNA\t{}".format(f, scaffold, position, chrom, output_cm, last_pattern[f]), file=crossfile)
    if crosstype is "Backcross":
        return last_pattern, cm
    return last_pattern

for chrom in range(1,22):
    bc_last_pattern = defaultdict(str)
    f1_last_pattern = defaultdict(str)
    all_last_pattern = defaultdict(str)

    bc_cm = defaultdict(float)
    
    with open('hybrids.map'+str(chrom)+'.txt') as m:
        for line in m:
            scaffold, position, map_cm, *gts = line.rstrip().split('\t')
            bc_this_pattern = defaultdict(str)
            f1_this_pattern = defaultdict(str)
            all_this_pattern = defaultdict(str)
            
            for i in range(len(gts)):
                bc_this_pattern[bc_header[i]] += gts[i]
                f1_this_pattern[hybrids[bc_header[i]]] += gts[i]
                all_this_pattern['Hybrid'] += gts[i]
            bc_last_pattern, bc_cm = output_pattern("Backcross", bc_this_pattern, bc_last_pattern, bc_file, bc_cm)
            f1_last_pattern  = output_pattern("F1", f1_this_pattern, f1_last_pattern, f1_file)
            all_last_pattern = output_pattern("All", all_this_pattern, all_last_pattern, all_file, map_cm)