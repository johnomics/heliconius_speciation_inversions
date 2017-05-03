#!/usr/bin/env python3

import argparse
from collections import defaultdict

chromosome_lengths = [0, 17206585, 9045316, 10541528, 9662098, 9908586, 14054175, 14308859, 9320449, 8708747, 17965481, 11759272, 16327298, 18127314, 9174305, 10235750, 10083215, 14773299, 16803890, 16399344, 14871695, 13359691]

parser=argparse.ArgumentParser(description='''
    -r recombinations
    -o output
    ''')

parser.add_argument('-r', '--recombinations', type=str, required=True)
parser.add_argument('-o', '--outputstub', type=str, required=True)
args = parser.parse_args()

speciesnames = ['Cydno','Hybrid','Melpomene']
crosses = {'CYDA':'Cydno', 'CYDB':'Cydno', 'CYDE':'Cydno', 'MEL1':'Melpomene', 'MEL3':'Melpomene', 'MEL4':'Melpomene', 'C6':'Hybrid', 'C26':'Hybrid', 'C29':'Hybrid', 'C1':'Hybrid', 'C94':'Hybrid'}
offspring = {'CYDA':95, 'CYDB':77, 'CYDE':125, 'MEL1':111, 'MEL3':122, 'MEL4':102, 'Hybrid':331,'Cydno':297,'Melpomene':335, 'C6':170, 'C26':88, 'C29':62, 'C1':5, 'C94':6}

recombinations = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
with open(args.recombinations) as recfile:
    for i, line in enumerate(recfile):
        if i == 0:
            hf = line.rstrip().split('\t')
            continue
        f = dict(zip(hf, line.rstrip().split('\t')))
        species, cross, chromosome, start, end = f['Species'], f['Cross'], int(f['Chromosome']), int(f['Start']), int(f['End'])
        recombinations[species][chromosome][(start, end)] += 1
        recombinations[cross][chromosome][(start, end)] += 1

breakpoints = defaultdict(lambda: defaultdict(list))
for name in recombinations:
    for chromosome in recombinations[name]:
        for start, end in recombinations[name][chromosome]:
            breakpoints[name][chromosome].append(start)
            breakpoints[name][chromosome].append(end)
        breakpoints[name][chromosome] = list(set(breakpoints[name][chromosome]))
        breakpoints[name][chromosome].sort()


def rec_count(win_start, win_end, rec_start, rec_end):
    rec_frac = 0
    rec_length = rec_end - rec_start + 1

    frac_start = max(rec_start, win_start)
    frac_end   = min(rec_end,   win_end  )

    if frac_start < frac_end:
        rec_frac = (frac_end - frac_start + 1) / rec_length

    return(rec_frac)

def write_cm_values(name, key, recombinations, breakpoints, outfile):

    for chromosome in sorted(breakpoints[name]):
        
        cm = 0
        print('{}\t{}\t0\t0\t0.000\t0.000'.format(key, chromosome), file=outfile)
        print('{}\t{}\t{}\t{}\t0.000\t0.000'.format(key, chromosome, 1, breakpoints[name][chromosome][0]), file=outfile)
    
        for i in range(len(breakpoints[name][chromosome])-1):
            win_start = breakpoints[name][chromosome][i]
            win_end   = breakpoints[name][chromosome][i+1]
            recs = 0
            for rec in recombinations[name][chromosome]:
                if win_start < rec[1] and win_end > rec[0]:
                    recfrac = rec_count(win_start, win_end, rec[0], rec[1])
                    recs += recfrac * recombinations[name][chromosome][rec]
            cmfrac = recs / offspring[name] * 100
            cm += cmfrac
            print('{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}'.format(key, chromosome, win_start, win_end, cmfrac, cm), file=outfile)
    
        print('{}\t{}\t{}\t{}\t0.000\t{:.3f}'.format(key, chromosome, win_end, chromosome_lengths[chromosome], cm), file=outfile)


speciesfile = open(args.outputstub + '.species.cm.tsv', 'w')
print("Species\tChromosome\tStart\tEnd\tcMFraction\tcM", file=speciesfile)
for species in speciesnames:
    write_cm_values(species, species, recombinations, breakpoints, speciesfile)
speciesfile.close()

crossesfile = open(args.outputstub + '.crosses.cm.tsv', 'w')
print("Species\tCross\tChromosome\tStart\tEnd\tcMFraction\tcM", file=crossesfile)
for cross in sorted(crosses):
    key = '{}\t{}'.format(crosses[cross], cross)
    write_cm_values(cross, key, recombinations, breakpoints, crossesfile)
crossesfile.close()


