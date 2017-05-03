#!/usr/bin/env python3

import argparse, sys, glob, os
from intervaltree import Interval, IntervalTree
from copy import deepcopy
import multiprocessing as mp
from termcolor import colored
from collections import defaultdict, namedtuple

sample_names = {'Cydno_maternal':'CF', 'Cydno_paternal':'CM', 'Melpomene_maternal':'MF', 'Melpomene_paternal':'MM'}
pacbio_names = {'Cydno_females':'CF', 'Cydno_males':'CM', 'Melpomene_females':'MF', 'Melpomene_males':'MM'}

chromosome_lengths = [0, 17206585, 9045316, 10541528, 9662098, 9908586, 14054175, 14308859, 9320449, 8708747, 17965481, 11759272, 16327298, 18127314, 9174305, 10235750, 10083215, 14773299, 16803890, 16399344, 14871695, 13359691]

class Inversion:
    def __init__(self, f):
        self.sample = f['Sample']
        self.species, self.sex = f['Sample'].split('_')
        self.species = self.species.title()
        self.chromosome = int(f['Chromosome1'])
        self.start = int(f['ChromPosition1'])
        self.end = int(f['ChromPosition2'])
        self.scaffold = f['Scaffold1']
        self.reads = int(f['Reads'])
        self.length = int(f['Length'])
        self.remainseq = int(f['RemainSeq'])
        self.initiallengths = f['InitialLengths'] if 'InitialLengths' in f else ''
        self.taillengths = f['TailLengths'] if 'TailLengths' in f else ''
        self.gaps = set()
        self.recs = set()
        self.hits = IntervalTree()
        self._status = ""
        
    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(self.species, self.sex, self.chromosome, self.start, self.end, self.scaffold, self.reads)
    
    @property
    def status(self):
        if self._status:
            return self._status

        status = ""
        for rec in self.recs:
            if rec.data[1] == self.species and rec.begin >= self.start and rec.end <= self.end:
                status = "Contains recombination"
        
        if not status:
            same, other, output = get_spanning_hits(IntervalTree(hits[self.chromosome].search(self.start, self.end)), Interval(self.start, self.end, self))
            if same:
                status += "Spanning hit"
        
        if not status:
            if self.length < args.ldminsize:
                status += "Shorter than LD threshold"

        if not status:
            status = "OK"
        
        self._status = status
        return status

class Hit:
    def __init__(self, samplename, tstart, tend, hstart, hend, thitlen, hhitlen, pctid, tscflen, hscflen, tcov, hcov, tdir, hdir, tname, hname):
        self.sample = samplename
        self.species, self.sex = self.sample.split('_')
        self.sex = "females" if self.sex == "maternal" else "males"
        self.chromosome = int(hname[5:7])
        self.tstart = tstart
        self.tend = tend
        
        self.hstart, self.hend = sorted([hstart, hend])
        
        self.thitlen = thitlen
        self.hhitlen = hhitlen
        self.pctid = pctid
        self.tscflen = tscflen
        self.hscflen = hscflen
        self.tcov = tcov
        self.hcov = hcov
        
        self.tdir = tdir
        self.hdir = hdir
        self.dir = tdir * hdir
        
        self.tname = tname
        self.hname = hname


class OutputLine:
    def __init__(self, groupid, invid, groupinvid, hittype, species, sex, chromosome, start, end, length, reads, pctid, dir, scaffold):
        self.groupid = groupid
        self.invid = invid
        self.groupinvid = groupinvid
        self.hittype = hittype
        self.species = species
        self.sex = sex
        self.chromosome = chromosome
        self.start = start
        self.end = end
        self.length = length
        self.reads = reads
        self.pctid = pctid
        self.dir = dir
        self.label = scaffold
    
    def __eq__(self, other):
        return self.hittype == other.hittype and self.species == other.species and self.sex == other.sex and self.chromosome == other.chromosome and self.start == other.start and self.end == other.end

    def __repr__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNA'.format(self.groupid, self.invid, self.groupinvid, self.hittype, self.species, self.sex, self.chromosome, self.start, self.end, self.length, self.reads, self.pctid, self.dir, self.label)

def spanning(tree, begin, end):
    return tree.search(begin).intersection(tree.search(end))

def overlap_edges(tree, begin, end):
    return tree.search(begin, end) - tree.search(begin, end, strict=True)

def get_map_regions(mapfilename):
    map_trees   = defaultdict(lambda: defaultdict(IntervalTree))

    with open(mapfilename,'r') as mapfile:
        for line in mapfile:
            if line.startswith("Species"):
                continue
    
            species, chromosome, start, end, cmfraction, cm = line.rstrip().split('\t')
            if species == "Hybrid":
                continue
            chromosome, start, end = int(chromosome), int(start), int(end)
            cmfraction, cm = float(cmfraction), float(cm)
            if cmfraction == 0:
                continue
            map_trees[species][chromosome][start:end] = (cm, species)
    
    return map_trees

def make_gap_trees(mapfilename):
    map_trees = get_map_regions(mapfilename)

    gap_trees = defaultdict(IntervalTree)
    for species in sorted(map_trees):
        for chromosome in sorted(map_trees[species]):
            cmlist = sorted(map_trees[species][chromosome])

            gap_trees[chromosome][1:cmlist[0].end] = species # First gap
            
            for i in range(len(cmlist)-1):
                gapstart  = cmlist[i].begin
                gapend    = cmlist[i+1].end
                gaplength = gapend - gapstart + 1
                gap_trees[chromosome][gapstart:gapend] = species

            gap_trees[chromosome][(cmlist[-1].begin-1):chromosome_lengths[chromosome]] = species # Last gap
 
            gap_trees[chromosome][(cmlist[-1][0]-1):chromosome_lengths[chromosome]] = species   # Last gap

    return gap_trees, map_trees

def write_gaps(gap_trees, outputstub):
    with open(outputstub + ".gaps.tsv", 'w') as outfile:
        print("Species\tChromosome\tGapStart\tGapEnd\tGapLength", file=outfile)
        for chromosome in gap_trees:
            for iv in sorted(gap_trees[chromosome]):
                print("{}\t{}\t{}\t{}\t{}".format(iv.data, chromosome, iv.begin, iv.end, iv.end-iv.begin), file=outfile)


def make_inversion_trees(inversionsfilename):
    total_inversions = valid_inversions = 0
    inv_trees = defaultdict(IntervalTree)
    samples = defaultdict(int)
    for i, line in enumerate(open(inversionsfilename)):
        if i == 0:
            hf = line.rstrip().split('\t')
            continue
        total_inversions += 1
        f = dict(zip(hf, line.rstrip().split('\t')))
        chr1, chr2, cp1, cp2 = int(f['Chromosome1']), int(f['Chromosome2']), int(f['ChromPosition1']), int(f['ChromPosition2'])
        if chr1 != chr2 or chr1 == 0 or f['Scaffold1'] != f['Scaffold2'] or cp2 <= cp1:
            continue
        valid_inversions += 1
        samples[f['Sample']] += 1
        inv_trees[chr1][cp1:cp2] = Inversion(f)
    
    return inv_trees

def make_candidate_groups(inv_trees, gap_trees, tentative, map_trees, hits):

    cand_groups = defaultdict(IntervalTree)
    ok_trees = defaultdict(IntervalTree)
    reject_trees = defaultdict(IntervalTree)
    statuses = defaultdict(lambda: defaultdict(int))
    
    for chromosome in sorted(inv_trees):
        for candidate in sorted(inv_trees[chromosome]):
            candidate.data.gaps = gap_trees[chromosome].search(candidate.begin, candidate.end)

            map_overlaps = set()
            for species in map_trees:
                map_overlaps |= map_trees[species][chromosome].search(candidate.begin, candidate.end)
            candidate.data.recs = map_overlaps

            if candidate.data.status == "OK":
                ok_trees[chromosome].add(candidate)
            else:
                reject_trees[chromosome].add(candidate)
                statuses[candidate.data.status][candidate.data.sample] += 1

        ok_tentative = IntervalTree()
        for candidate in ok_trees[chromosome]:
            tenpclen = candidate.length() * 0.1
            for ct in tentative[chromosome].search(candidate.begin, candidate.end):
                reject = False
                for ctok in ok_trees[chromosome].search(candidate.begin, candidate.end):
                    if ct.begin == ctok.begin and ct.end == ctok.end and ct.data.species == ctok.data.species and ct.data.sex == ctok.data.sex:
                        reject = True
                        continue
                if not reject and \
                   abs(ct.begin-candidate.begin) < tenpclen and abs(ct.end-candidate.end) < tenpclen and \
                   ct.data.status == "OK":
                    if ct.data.status == "OK":
                        ct.data._status = "Tentative"
                        ok_tentative.add(ct)

        for okt in ok_tentative:
            ok_trees[chromosome].add(okt)
        
        for ok in ok_trees[chromosome]:
            statuses[ok.data.status][ok.data.sample] += 1

        cand_groups[chromosome] = ok_trees[chromosome].copy()
        cand_groups[chromosome].merge_overlaps()

    return cand_groups, ok_trees, reject_trees

def load_transitions(transitionfile):
    hmel2 = defaultdict(IntervalTree)
    hmel2o = defaultdict(int)
    chrom_to_hmel2 = defaultdict(IntervalTree)
    
    for line in open(transitionfile):
        if line.startswith("Chromosome"):
            continue
        f = line.rstrip().split('\t')
        chromosome, chromstart, chromend, pool,       pooltype, poolid,    hmel2scaffold, hmel2start, hmel2end,  orientation, hmel2refined, length,     hmel2ordered, hmel2ostart, hmel2oend = \
        int(f[0]),  int(f[1]),  int(f[2]), int(f[3]), f[4],     int(f[5]), f[6],          int(f[7]),  int(f[8]), f[9],        f[10],        int(f[11]), f[12],        int(f[13]),  int(f[14])

        if chromosome == 0:
            continue
        
        if hmel2ordered not in hmel2o:
            hmel2o[hmel2ordered] = chromstart
        
        hmel2[hmel2scaffold][hmel2start:hmel2end] = (chromosome, chromstart, chromend, orientation)
        chrom_to_hmel2[chromosome][chromstart:chromend] = (hmel2scaffold, hmel2start, hmel2end, orientation)

    return hmel2, hmel2o, chrom_to_hmel2

def load_trio_alignments(args, hmel2o):
    
    hits = defaultdict(IntervalTree)
    coordsnames = glob.glob(args.coordsglob)
    for coordsname in coordsnames:
        print("Loading " + coordsname)
        samplename = coordsname.split(os.sep)[-1].split('.')[0]
        for line in open(coordsname):
            if line.endswith('o\n'):
                f = line.rstrip().split('\t')
                tstart,  tend,    hstart, hend = int(f[0]),  int(f[1]),  int(f[2]),   int(f[3])
                thitlen, hhitlen, pctid        = int(f[4]),  int(f[5]),  float(f[6])
                tscflen, hscflen, tcov,   hcov = int(f[7]),  int(f[8]),  float(f[9]), float(f[10])
                tdir,    hdir,    tname, hname = int(f[11]), int(f[12]), f[13],       f[14]

                if thitlen < args.minhitlength:
                    continue
                
                tname = tname[1:] # Remove erroneous leading > from FASTA name
                locname = samplename + '_' + tname
                chromosome = int(hname[5:7])
                hstart     = hstart + hmel2o[hname] - 1
                hend       = hend   + hmel2o[hname] - 1
                hstart, hend = sorted([hstart, hend])
                hits[chromosome].addi(hstart, hend, Hit(samplename, tstart, tend, hstart, hend, thitlen, hhitlen, pctid,
                                         tscflen, hscflen, tcov, hcov, tdir, hdir, tname, hname))

    return hits



def get_gff_chrom_position(scaffold, start, end, hmel2):
    hmel2_overlaps = hmel2[scaffold].search(start, end)
    if len(hmel2_overlaps) != 1:
        return None, None, None

    hmel2ov = list(hmel2_overlaps)[0]
    hmel2scfstart = hmel2ov.begin
    hmel2chrom, hmel2chromstart, hmel2chromend, hmel2orient = hmel2ov.data
    if hmel2orient == '+':
        chromstart = start - hmel2scfstart + hmel2chromstart - 1
        chromend   = end   - hmel2scfstart + hmel2chromstart - 1
    else:
        chromstart = hmel2chromend - (start - hmel2scfstart)
        chromend   = hmel2chromend - (end - hmel2scfstart)
    chromstart, chromend = sorted([chromstart, chromend])
    
    return hmel2chrom, chromstart, chromend

def get_hmel2_position(chromosome, start, end, chrom_to_hmel2):

    hmel2_overlaps = chrom_to_hmel2[chromosome].search(start, end)
    if len(hmel2_overlaps) != 1:
        return None, None, None

    hmel2ov = list(hmel2_overlaps)[0]
    chromstart = hmel2ov.begin
    hmel2scaffold, hmel2scfstart, hmel2scfend, hmel2orient = hmel2ov.data
    if hmel2orient == '+':
        hmel2start = start - chromstart + hmel2scfstart - 1
        hmel2end   = end   - chromstart + hmel2scfstart - 1
    else:
        hmel2start = hmel2scfend - (start - chromstart)
        hmel2end   = hmel2scfend - (end - chromstart)

    hmel2start, hmel2end = sorted([hmel2start, hmel2end])
    
    return hmel2scaffold, hmel2start, hmel2end


def load_gff(gfffile, hmel2=None):

    gff_trees = defaultdict(IntervalTree)
    for line in open(gfffile):
        if line.startswith('#'):
            continue
        scaffold, origin, featuretype, start, end, score, strand, frame, attributes = line.rstrip().split('\t')
        chromosome = int(scaffold[5:7])
        if chromosome == 0:
            continue

        gffchrom, gffstart, gffend = chromosome, int(start), int(end)
        if hmel2 is not None:
            gffchrom, gffstart, gffend = get_gff_chrom_position(scaffold, gffstart, gffend, hmel2)
        
        if gffchrom is None:
            continue
        
        gff_trees[gffchrom][gffstart:gffend] = (featuretype, attributes)
        
    return gff_trees


def load_agp(agpfile):
    agp_trees = defaultdict(IntervalTree)
    for line in open(agpfile):
        f = line.rstrip().split('\t')
        chromosome, start, end = int(f[0]), int(f[1]), int(f[2])
        data = f[5] if f[4] == 'W' else 'Gap'
        agp_trees[chromosome][start:end+1] = data
    
    return(agp_trees)


def get_spanning_hits(hits, inv_iv, tsv_hits=None, group_details=None):
    spanningtree = spanning(hits, inv_iv.begin, inv_iv.end)
    sample_hits = defaultdict(list)
    
    for hit in spanningtree:
        leftlen = inv_iv.begin - hit.data.hstart + 1
        rightlen   = hit.data.hend - inv_iv.end + 1
        leftlenratio = leftlen / inv_iv.length()
        rightlenratio   = rightlen   / inv_iv.length()
        
        if leftlenratio < 0.5 or rightlenratio < 0.5:
            continue

        sample_hits[hit.data.sample].append((leftlen, rightlen, leftlenratio, rightlenratio))

        if tsv_hits and group_details:
            tsv_hits.append(OutputLine(group_details[0], group_details[1], group_details[2], "Spanning", hit.data.species, hit.data.sex, group_details[3], hit.data.hstart, hit.data.hend, hit.data.hhitlen, 'NA', hit.data.pctid, hit.data.dir, hit.data.tname))
    
    same = other = 0
    for s in sample_hits:
        if len(sample_hits[s]) == 0:
            continue
        if inv_iv.data.species in s:
            same += 1
        else:
            other += 1
    return same, other

def add_edge_hits(edgetree, sample_hits, point, length, inverted_scaffolds):
    for hit in edgetree:
        leftlen    = point - hit.data.hstart + 1
        rightlen   = hit.data.hend - point + 1
        leftlenratio  = leftlen  / length
        rightlenratio = rightlen / length
        if hit.data.tname not in inverted_scaffolds[hit.data.sample]:
            sample_hits[hit.data.sample].append((leftlen, rightlen, leftlenratio, rightlenratio))

def write_edge_tsv(edgetree, begin_tree, end_tree, tsv_hits, group_details, inverted_scaffolds):
    scaffolds = {}
    for iv in begin_tree:
        scaffolds[iv.data.tname] = 1
    for iv in end_tree:
        scaffolds[iv.data.tname] = 1


    for hit in edgetree:
        if hit.data.tname in scaffolds and hit.data.tname not in inverted_scaffolds[hit.data.sample]:
            tsv_hits.append(OutputLine(group_details[0], group_details[1], group_details[2], "Edge", hit.data.species, hit.data.sex, group_details[3], hit.data.hstart, hit.data.hend, hit.data.hhitlen, 'NA', hit.data.pctid, hit.data.dir, hit.data.tname))
    

def get_edge_hits(hits, inv_iv, tsv_hits, group_details, inverted_scaffolds):

    begintree = hits.search(inv_iv.begin)
    endtree   = hits.search(inv_iv.end)
    begin_edgetree = begintree - endtree
    end_edgetree = endtree - begintree

    sample_hits = defaultdict(list)

    add_edge_hits(begin_edgetree, sample_hits, inv_iv.begin, inv_iv.length(), inverted_scaffolds)
    add_edge_hits(end_edgetree, sample_hits, inv_iv.end, inv_iv.length(), inverted_scaffolds)
    
    write_edge_tsv(hits.search(inv_iv.begin, inv_iv.end), begin_edgetree, end_edgetree, tsv_hits, group_details, inverted_scaffolds)
    
    return


def get_outer_hits(outer_hits, scaffold, sample, hitdir):
    scf_dir_tree = IntervalTree()
    for outer_hit in outer_hits:
        if scaffold == outer_hit.data.tname and sample == outer_hit.data.sample and hitdir != outer_hit.data.dir:
            scf_dir_tree.add(outer_hit)
    return scf_dir_tree

def get_inverted_hits(hits, inv_iv, tsv_hits, group_details):

    hits_by_scaffold=defaultdict(lambda: defaultdict(lambda: defaultdict(IntervalTree)))
    for hit in hits.search(inv_iv.begin, inv_iv.end):
        hits_by_scaffold[hit.data.tname][hit.data.sample][hit.data.dir].add(hit)
    
    left_outer_hits = hits.search(inv_iv.begin - args.extendregion, inv_iv.begin)
    right_outer_hits = hits.search(inv_iv.end, inv_iv.end + args.extendregion)

    sample_hits = defaultdict(list)

    inverted_scaffolds = defaultdict(lambda: defaultdict(int))
    for scaffold in sorted(hits_by_scaffold):
        for sample in sorted(hits_by_scaffold[scaffold]):
            for hitdir in sorted(hits_by_scaffold[scaffold][sample]):
                withinrange = hits_by_scaffold[scaffold][sample][hitdir].range()
                
                scf_dir_left = get_outer_hits(left_outer_hits, scaffold, sample, hitdir)
                scf_dir_right = get_outer_hits(right_outer_hits, scaffold, sample, hitdir)

                num_left, num_right = len(scf_dir_left), len(scf_dir_right)
                if num_left + num_right > 0:
                    sample_hits[sample].append((withinrange.begin, withinrange.end, len(scf_dir_left), len(scf_dir_right)))
                    
                    for hit in [h for tree in [hits_by_scaffold[scaffold][sample][hitdir], scf_dir_left, scf_dir_right] for h in tree]:
                        tsv_hits.append(OutputLine(group_details[0], group_details[1], group_details[2], "Inverted", hit.data.species, hit.data.sex, group_details[3], hit.data.hstart, hit.data.hend, hit.data.hhitlen, 'NA', hit.data.pctid, hit.data.dir, hit.data.tname))
                        inverted_scaffolds[sample][hit.data.tname] = 1
                        

    same = other = 0
    for s in sample_hits:
        if len(sample_hits[s]) == 0:
            continue
        if inv_iv.data.species in s:
            same += 1
        else:
            other += 1
    return same, other, inverted_scaffolds


def summarise_candidate(inv_iv, statuses, spanning_same, spanning_other, inverted_same, inverted_other):
    status = 'None'
    if spanning_same and not inverted_same and not inverted_other:
        status = 'Reject'
    if not spanning_same and inverted_same and not inverted_other:
        status = 'Accept'
    if not spanning_same and not spanning_other and inverted_same and inverted_other:
        status = 'Misassembly'
    statuses[status][pacbio_names[inv_iv.data.species + '_' + inv_iv.data.sex]] += 1
    return status


def get_group_samples(samples, group_num, inv_num, in_group_num, chromosome, candi, inv_iv):

    for sn in sorted(sample_names):
        shortsn = sample_names[sn]
        if shortsn == pacbio_names[inv_iv.data.species + '_' + inv_iv.data.sex]:
            samples += shortsn + '_'
    return samples

def get_status(group_samples, trio_statuses):
    pbhoney_status = 'Single'
    if ('CF' in group_samples or 'CM' in group_samples) and ('MF' in group_samples or 'MM' in group_samples):
        pbhoney_status = 'Misassembly'
    elif 'CM' in group_samples and 'CF' in group_samples:
        pbhoney_status = 'Cydno'
    elif 'MM' in group_samples and 'MF' in group_samples:
        pbhoney_status = 'Melpomene'

    status = '????'
    if pbhoney_status == 'Single':
        trio_status = list(trio_statuses)[0]
        if trio_status == 'Misassembly':
            status = 'Misassembly_Trio'
        elif trio_status == 'None':
            if 'CF' in group_samples or 'CM' in group_samples:
                status = 'Cydno_Single'
            else:
                status = 'Melpomene_Single'
        elif trio_status == 'Accept':
            if 'CF' in group_samples or 'CM' in group_samples:
                status = 'Cydno_Both'
            elif 'MF' in group_samples or 'MM' in group_samples:
                status = 'Melpomene_Both'
            else:
                status = '????'
        else:
            status = '????'
    elif pbhoney_status == 'Misassembly':
        if 'Misassembly' in trio_statuses:
            status = 'Misassembly_Both'
        else:
            status = 'Misassembly_PBHoney'
    else: # pbhoney status is Cydno or Melpomene
        if 'Misassembly' in trio_statuses:
            status = 'Misassembly_Trio'
        elif 'Accept' in trio_statuses:
            status = pbhoney_status + '_Both'
        elif len(trio_statuses) == 1 and 'None' in trio_statuses:
            status = pbhoney_status + '_PBHoney'
        else:
            status = '????'
    return status


def write_tsv(group, group_num, group_status, chromosome, tsv_inversions, tsv_hits, plottsv, gff_trees, chrom_to_hmel2, agp_trees, repeats):
    
    print('{}\t0\t0\tGroup\tAll\t\t{}\t{}\t{}\t{}\tNA\tNA\tNA\tAll\t{}\t{}'.format(group_num, chromosome, group.begin, group.end, group.length(), group_status, group_status), file=plottsv)

    for inv in tsv_inversions:
        print(inv + "\t" + group_status, file=plottsv)
    
    prevhit = None
    for hit in sorted(tsv_hits, key = lambda x: (x.hittype, x.species, x.sex, x.chromosome, x.start)):
        if prevhit and hit != prevhit:
            print(repr(hit) + "\t" + group_status, file=plottsv)
        prevhit = hit
        
    for gff_feature in sorted(gff_trees[chromosome].search(group.begin, group.end)):
        print('{}\t0\t0\tFeature\tAll\t\t{}\t{}\t{}\t{}\tNA\tNA\tNA\t{}\t{}\t{}'.format(group_num, chromosome, gff_feature.begin, gff_feature.end, gff_feature.length(), gff_feature.data[0], gff_feature.data[1], group_status), file=plottsv)


    for agp_part in sorted(agp_trees[chromosome].search(group.begin-group.length()/2, group.end+group.length()/2)):
        print('{}\t0\t0\tHmel2Part\tAll\t\t{}\t{}\t{}\t{}\tNA\tNA\tNA\t{}\tNA\t{}'.format(group_num, chromosome, agp_part.begin, agp_part.end, agp_part.length(), agp_part.data, group_status), file=plottsv)
    
    for repeat in sorted(repeats[chromosome].search(group.begin-group.length()/2, group.end+group.length()/2)):
        print('{}\t0\t0\tRepeat\tAll\t\t{}\t{}\t{}\t{}\tNA\tNA\tNA\t{}\tNA\t{}'.format(group_num, chromosome, repeat.begin, repeat.end, repeat.length(), repeat.data, group_status), file=plottsv)


def get_args():
    parser=argparse.ArgumentParser(description='''Find and classify candidate inversions
        -m map
        -i inversions
        -c coordsglob
        -o output
        -s transitions
        -l minhitlength
        -t tentativeinversions
        -g gff
        -a agp
        -r repeats
        -d ldminsize
        ''')
    
    parser.add_argument('-m', '--map', type=str, required=True)
    parser.add_argument('-i', '--inversions', type=str)
    parser.add_argument('-o', '--outputstub', type=str, required=True)
    parser.add_argument('-c', '--coordsglob', type=str)
    parser.add_argument('-s', '--transitions', type=str)
    parser.add_argument('-l', '--minhitlength', type=int, default=1000)
    parser.add_argument('-e', '--extendregion', type=int, default=10000)
    parser.add_argument('-t', '--tentativeinversions', type=str)
    parser.add_argument('-g', '--gff', type=str)
    parser.add_argument('-a', '--agp', type=str)
    parser.add_argument('-r', '--repeats', type=str)
    parser.add_argument('-d', '--ldminsize', type=int, default=1000)
    
    return parser.parse_args()
    
args = get_args()

print("Making and writing gaps")

gap_trees, map_trees = make_gap_trees(args.map)
write_gaps(gap_trees, args.outputstub)

if not args.inversions:
    sys.exit()

print("Loading PBHoney inversions")
inv_trees = make_inversion_trees(args.inversions)

if args.tentativeinversions:
    print("Loading tentative inversions")
    tentative = make_inversion_trees(args.tentativeinversions)
else:
    tentative = defaultdict(IntervalTree)

if not args.transitions and (args.gff or args.coordsglob):
    print("Need Hmel2 transition file (-s)")
    sys.exit()

hmel2, hmel2o, chrom_to_hmel2 = load_transitions(args.transitions)

if args.gff:
    print("Loading transcriptome")
    gff_trees = load_gff(args.gff, hmel2)
else:
    gff_trees = defaultdict(IntervalTree)

if args.agp:
    print("Loading AGP file")
    agp_trees = load_agp(args.agp)
else:
    agp_trees = defaultdict(IntervalTree)


if args.repeats:
    print("Loading repeats")
    repeats = load_gff(args.repeats)
else:
    repeats = defaultdict(IntervalTree)

if not args.coordsglob:
    sys.exit()

print("Loading trio alignments")
hits = load_trio_alignments(args, hmel2o)

print("Making candidate groups")
cand_groups, ok_trees, reject_trees = make_candidate_groups(inv_trees, gap_trees, tentative, map_trees, hits)

group_num = inv_num = 0
summary = defaultdict(lambda: defaultdict(int))

plottsv = open(args.outputstub + ".details.tsv", 'w')
print("GroupID\tInvID\tGroupInvID\tHitType\tSpecies\tSex\tChromosome\tStart\tEnd\tLength\tReads\tPctId\tDir\tLabel\tStatus\tGroupStatus", file=plottsv)

for chromosome in sorted(cand_groups):
    for group in sorted(cand_groups[chromosome]):
        group_num += 1
        in_group_num = 0
        group_samples = ''
        tsv_hits = []
        tsv_inversions = []
        trio_statuses = defaultdict(lambda: defaultdict(int))

        group_inversions = sorted(ok_trees[chromosome].search(group.begin, group.end, strict=True), key=lambda x:(x.data.species, x.data.sex, x.begin))

        for inv_iv in group_inversions:
            inv_num += 1
            in_group_num += 1

            group_samples = get_group_samples(group_samples, group_num, inv_num, in_group_num, chromosome, group, inv_iv)
            group_details = (group_num, inv_num, in_group_num, chromosome)

            spanning_same, spanning_other                     = get_spanning_hits(hits[chromosome], inv_iv, tsv_hits, group_details)
            inverted_same, inverted_other, inverted_scaffolds = get_inverted_hits(hits[chromosome], inv_iv, tsv_hits, group_details)
            get_edge_hits(hits[chromosome], inv_iv, tsv_hits, group_details, inverted_scaffolds)

            inv_status = summarise_candidate(inv_iv, trio_statuses, spanning_same, spanning_other, inverted_same, inverted_other)
        
            tsv_inversions.append('{0}\t{1}\t{2}\tPBHoney\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{3}_{4} Candidate\t{12}'.format(group_num, inv_num, in_group_num, inv_iv.data.species,  inv_iv.data.sex, chromosome, inv_iv.begin, inv_iv.end, inv_iv.length(), inv_iv.data.reads, inv_iv.data.initiallengths, inv_iv.data.taillengths, inv_status))

        group_status = get_status(group_samples, trio_statuses)
        write_tsv(group, group_num, group_status, chromosome, tsv_inversions, tsv_hits, plottsv, gff_trees, chrom_to_hmel2, agp_trees, repeats)

        summary[group_status][group_samples[:-1]] += 1

    for reject in sorted(reject_trees[chromosome]):
        print('0\t0\t0\tReject\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNA\t{}\tReject'.format(reject.data.species, reject.data.sex, chromosome, reject.begin, reject.end, reject.length(), reject.data.reads, reject.data.initiallengths, reject.data.taillengths, reject.data.status), file=plottsv)

plottsv.close()