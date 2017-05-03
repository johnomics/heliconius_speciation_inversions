This repository documents analyses for the manuscript "No evidence for maintenance of a sympatric _Heliconius_ species barrier by chromosomal inversions". The bespoke scripts are presented for transparency and will require editing to be applied to other datasets. This Github repository contains only the code written for the manuscript. All data files can be found in the [Dryad repository for the manuscript](http://dx.doi.org/10.5061/dryad.kv920). Please contact John Davey (johnomics on Twitter and Gmail) for help.

Any file from the Hmel2 distribution (version 2 of the _H. melpomene_ genome assembly) refers to the distribution of [13 October 2015](http://butterflygenome.org/sites/default/files/Hmel2-0_Release_20151013.tgz), available from [butterflygenome.org](http://butterflygenome.org) and [LepBase](http://ensembl.lepbase.org/Heliconius_melpomene_hmel2/Info/Index), but the relevant files have also been included here in directory `Hmel2`.

# Linkage mapping


## Within-species maps

Within-species maps for _H. melpomene_ and _H. cydno_ were constructed using [Lep-MAP2](https://sourceforge.net/projects/lepmap2/) and modules from [Lep-MAP3](https://sourceforge.net/projects/lep-map3/) (see Methods for full details). Six maps were constructed, for melpomene crosses 1, 2, and 3 (MEL1, MEL3, and MEL4 respectively) and cydno crosses 1, 2, and 3 (CYDA, CYDB, and CYDE respectively) using the same procedure. MEL1 is used for the following example commands.

Lep-MAP outputs the following files for each cross:
```
MEL1.call.gz   - genotype posteriors
MEL1_snps.txt  - SNP positions
MEL1_id.err    - unique marker patterns
MEL1_id.txt    - marker patterns assigned to SNPs
```

And 21 map files called `MEL1.order*.txt` for chromosomes 1 to 21.

These maps were cleaned manually; `clean_Lep-MAP_output.pl` was run once to identify inconsistencies, which were then added to `MEL1.errors.tsv`, and then run again to clean the maps. The script also reverses chromosomes listed in the manually compiled file `MEL1.reverse.txt` to ensure linkage maps are ordered as per Hmel2 chromosomes. It produces output files `MEL1.markers.tsv`. `-c` specifies the number of chromosomes, `-s` the species name and `-p` the cross name.

```
scripts/clean_Lep-MAP_output.pl -c 21 -s Melpomene -p MEL1 -o
```

Cleaned maps were then processed by `assign_snps_to_markers.pl`, which reassigned to SNPs to improve coverage across the genome, and attempted to find intermediate markers between markers separated by multiple recombinations, using Lep-MAP3 module JoinSingles2. `-i` searches for intermediates, and `-j` specifies the path to JoinSingles2.

```
scripts/assign_snps_to_markers.pl -m MEL1.markers.tsv -p MEL1.call.gz -s Melpomene -c MEL1 -i -j ~/bin/Lep-MAP3
```

This script produces output files `MEL1.intermediates.markers.updated.tsv` and `MEL1.intermediates.map.tsv`. Segregation patterns were transferred from the `markers.updated.tsv` file to `map.tsv` file with `add_patterns_to_snp_file.py`:

```
scripts/add_patterns_to_snp_file.py -m MEL1.intermediates.markers.updated.tsv -s MEL1.intermediates.map.tsv -o MEL1.intermediates.map.patterns.tsv
```

All six `intermediates.map.patterns.tsv` were then compiled into one map, stripping out all headers except the first:

```
cat *.patterns.tsv | awk 'NR==1 || $1 !~ "Species"' > linkage_map.melpomene_cydno.Hmel2.tsv
```


## Hybrid maps

Hybrid _H. cydno_ x _H. melpomene_ maps were constructed using Lep-MAP3 from ParentCallGrandparents output `hybrids.call.gz` (see Methods for full details). The resulting maps are `linkage_maps/hybrids/hybrids.map*.txt` where `*` is chromosomes 1 to 21 and `linkage_maps/hybrids/header.txt` lists the crosses and individuals.

The complete map was separated by cross by `split_hybrids_by_cross.py`, which loads `header.txt` and each chromosome map in turn and produces three files containing markers grouped by families: `linkage_map.hybrids.backcross.Hmel2.tsv` contains the markers split by the 18 backcross families, `linkage_map.hybrids.f1s.Hmel2.tsv` contains the markers split by 4 F1 fathers, and `linkage_map.hybrids.all.Hmel2.tsv` contains the markers with all families merged together to produce one hybrid map.
```
scripts/split_hybrids_by_cross.py
```


# Genome ordering

[Hmel2](http://www.g3journal.org/content/6/3/695.long) was edited and ordered based on manual inspection of the new _H. melpomene_ linkage maps, with the new scaffold order recorded in manually created file `Hmel2.linkage_map_order.tsv`. In this file, scaffolds are grouped into pools and ordered along chromosomes. Pools can contain one scaffold with known orientation ('ok'), one scaffold with unknown orientation ('orient') or many scaffolds with unknown order and orientations ('order'). A small number of scaffolds were split (eg Hmel203009 became Hmel203009a and Hmel203009b) and some scaffolds were moved to different chromosomes (eg Hmel208049, previously on chromosome 8, is now Hmel206_08049 on chromosome 6).

Scaffolds were ordered further by aligning PacBio reads to Hmel2 and identifying reads overlapping scaffold ends. Scaffold connections were summarised with script `find_pacbio_scaffold_overlaps.py`, which searches in the directories specified by `-p` for sorted BAM files named, for example, `melpomene_males/melpomene_males.sorted.bam`:

```
scripts/find_pacbio_scaffold_overlaps.py -p melpomene_males,melpomene_females -s genome/Hmel2.linkage_map_order.tsv -a > genome/Hmel2.bridges.all.txt
```

Scaffolds were then ordered by manual inspection of `Hmel2.bridges.all.txt` to make file `Hmel2.pacbio_order.tsv`.

This new order was then used to generate a set of updated genomes and genome files from the Hmel2 assembly `Hmel2.fa` with script `reorder_Hmel2.py`:

```
scripts/reorder_Hmel2.py -s Hmel2.pacbio_order.tsv -f Hmel2/Hmel2.fa -o Hmel2
```

This script generates the following files:
```
Hmel2_refined.fa             - Hmel2 scaffolds split and renamed
Hmel2_ordered.fa             - super-scaffolded Hmel2 scaffolds as reported in paper (38 scaffolds placed on chromosomes)
Hmel2.transitions.tsv        - master list of chromosome positions for Hmel2, Hmel2_refined and Hmel2_ordered scaffold pieces
Hmel2_refined.chromosome.agp - chromosome positions for Hmel2_refined scaffolds
Hmel2_ordered.chromosome.agp - chromosome positions for Hmel2_ordered scaffolds
```

The script `transfer_Hmel2_agp.py` reports chromosome positions for Hmel2 contigs as per the new scaffold order.

```
scripts/transfer_Hmel2_agp.py -t Hmel2.transitions.tsv -a Hmel2/Hmel2_scaffolds.agp > Hmel2.chromosome.agp
```

Chromosome positions for Hmel1.1 scaffolds were calculated with script `transfer_Hmel1_positions.py`:
```
scripts/transfer_Hmel1_positions.py -i Hmel2/Hmel2_transfer_new.tsv -t genome/Hmel2.transitions.tsv > genome/Hmel1_chromosome_positions.tsv
```


# Recombination analysis

Within-species maps were transferred from Hmel2 positions to chromosome positions:
```
scripts/transfer_Hmel2_positions.py -t genome/Hmel2.transitions.tsv -i linkage_maps/linkage_map.melpomene_cydno.Hmel2.tsv > linkage_maps/linkage_map.melpomene_cydno.tsv
```

Hybrid maps were then transferred to chromosome positions:

```
scripts/transfer_Hmel2_positions.py -t genome/Hmel2.transitions.tsv -i linkage_maps/hybrids/linkage_map.hybrids.backcross.Hmel2.tsv -r > linkage_maps/hybrids/linkage_map.backcross.hybrids.tsv
scripts/transfer_Hmel2_positions.py -t genome/Hmel2.transitions.tsv -i linkage_maps/hybrids/linkage_map.hybrids.f1s.Hmel2.tsv -r > linkage_maps/hybrids/linkage_map.f1s.hybrids.tsv
scripts/transfer_Hmel2_positions.py -t genome/Hmel2.transitions.tsv -i linkage_maps/hybrids/linkage_map.hybrids.all.Hmel2.tsv -r > linkage_maps/hybrids/linkage_map.all.hybrids.tsv
```


Within-species and hybrid linkage maps were compiled together with `cat` and multiple header lines removed with `vi`:

```
cat linkage_map.melpomene_cydno.tsv hybrids/linkage_map.backcross.hybrids.tsv > linkage_map.tsv
vi linkage_map.tsv
```

Recombinations for each offspring were compiled with `collapse_maps.py`, filtering a manually compiled list of map inconsistencies in `linkage_map.errors.tsv`, and outputting `linkage_map.recombinations.tsv` and `linkage_map.crosses.map.tsv`:

```
scripts/collapse_maps.py -m linkage_maps/linkage_map.tsv -e linkage_maps/linkage_map.errors.tsv -o linkage_map
```

Permutation tests for recombination rate differences were run by script `permute_recombination_differences.R`, loading the output of `collapse_maps.py` and running on 1 Mb windows (`-w 1000000`) in 100 kb steps (`-s 100000`) over 270,000 permutations (`-p 270000`) and 50 threads (`-t 50`):

```
scripts/permute_recombination_differences.R -r linkage_map.recombinations.tsv -o linkage_map.permutations -w 1000000 -s 100000 -p 270000 -t 50
```

This script generates p-values for differences between the species for each window in `linkage_map.permutations.cm_differences.tsv`.


## Figure S2

Plots of genetic and physical maps in `FigureS2.pdf` were created with `plot_species_crosses.R`:
```
scripts/plot_species_crosses.R -t genome/Hmel2.transitions.tsv -s linkage_maps/linkage_map.crosses.map.tsv -m linkage_maps/linkage_map.tsv -a linkage_maps/hybrids/linkage_map.all.hybrids.tsv -o linkage_maps.plot
```

The resulting 21 PDFs for each chromosome were merged together with [PDFsam](http://www.pdfsam.org).


## Figure 2 and Figure S4

CentiMorgan values for species and crosses were calculated with `calculate_cm_values.py`, producing `linkage_map.species.cm.tsv` and `linkage_map.crosses.cm.tsv`:

```
scripts/calculate_cm_values.py -r linkage_map.recombinations.tsv  -o linkage_map
```

The `plot_marey_maps.R` script loads `linkage_map.species.cm.tsv` and `linkage_map.crosses.cm.tsv` and outputs `linkage_map.crosses.map.tsv` (used for Figures S11-S17, see below), `linkage_map.species.marey.pdf` (Figure 2) and `linkage_map.crosses.marey.pdf` (Figure S4).

```
scripts/plot_marey_maps.R -i linkage_map
```


## Figure S5

The `calculate_recombination_rates.R` script produces `linkage_map.recombination_rate.species.pdf`, which is `FigureS5.pdf`, and `linkage_map.recombination_rate.windows.species.tsv`. This command calculates rates for 1 Mb windows (`-w 1000000`) in 100 kb steps (`-s 100000`) for 10,000 iterations (`-i 10000`) using 50 threads (`-t 50`).

```
scripts/calculate_recombination_rates.R -r linkage_map.recombinations.tsv -o linkage_map -w 1000000 -s 100000 -i 10000 -t 50
```

## Figure S3

GC content and number of restriction sites were calculated with `count_gc_radsites_per_window.py`, using 1 Mb windows (`-w 1000000`) in 100 kb steps (`-s 100000`) for the PstI recognition site (`-r CTGCAG`):
```
scripts/count_gc_radsites_per_window.py -g genome/Hmel2_ordered.fa -w 1000000 -s 100000 -r CTGCAG > genome/Hmel2_ordered.gc.psti.tsv
```

The `calculate_snp_densities.R` script produces `linkage_map.snp_density.species.pdf`, which is `FigureS3.pdf`, and `linkage_map.snp_density.windows.species.tsv`. It also calculates correlations between SNP density, PstI sites, recombination rate and GC content.

```
scripts/calculate_snp_densities.R -r linkage_maps/linkage_map.crosses.map.tsv -o linkage_map -f genome/Hmel2_ordered.gc.psti.tsv -c linkage_maps/linkage_map.recombination_rate.windows.species.tsv
```


# Inversions

## PBHoney

Raw PBHoney output for each sample is contained in the following files:
```
cydno_females.PBHoney.tails
cydno_males.PBHoney.tails
melpomene_females.PBHoney.tails
melpomene_males.PBHoney.tails
cydno_females.PBHoney.tentative.tails
cydno_males.PBHoney.tentative.tails
melpomene_females.PBHoney.tentative.tails
melpomene_males.PBHoney.tentative.tails
```

These were merged together with the script `compile_tails.py`:

```
scripts/compile_tails.py -t genome/Hmel2.transitions.tsv -i PBHoney.tails
scripts/compile_tails.py -t genome/Hmel2.transitions.tsv -i PBHoney.tentative.tails
```

These runs produce output files `PBHoney.tails.inversions.tsv` and `PBHoney.tentative.tails.inversions.tsv`.

## Trio assemblies

Trio assemblies are in the following files:
```
inversions/Cydno_maternal.trio_assembly.fa
inversions/Cydno_paternal.trio_assembly.fa
inversions/Melpomene_maternal.trio_assembly.fa
inversions/Melpomene_paternal.trio_assembly.fa
```

Nucmer alignments of trio assemblies to Hmel2_ordered are in the following files:
```
inversions/Cydno_maternal.Hmel2_ordered.nucmer.coords
inversions/Cydno_paternal.Hmel2_ordered.nucmer.coords
inversions/Melpomene_paternal.Hmel2_ordered.nucmer.coords
inversions/Melpomene_maternal.Hmel2_ordered.nucmer.coords
```

## Analysis

`process_inversions.py` identifies gaps in the linkage map, producing output file `inversions.gaps.tsv`, and then identifies, filters and summarises the PBHoney candidate inversions using trio assembly, repeat and genome information. `Hmel2.gff` is the annotation in the Hmel2 distribution. `Hmel2_ordered.repeats.gff` contains the positions of repeats in Hmel2 identified by RepeatMasker using `-species Heliconius`. All output is written to `inversions.details.tsv`.

```
scripts/process_inversions.py -m linkage_maps/linkage_map.species.cm.tsv -i inversions/PBHoney.tails.inversions.tsv -c 'inversions/*coords' -s genome/Hmel2.transitions.tsv -g Hmel2/Hmel2.gff -a genome/Hmel2.chromosome.agp -r genome/Hmel2_ordered.repeats.gff -t inversions/PBHoney.tentative.tails.inversions.tsv -o inversions
```

## Figure 3


Figure 3 was created with script `Figure3.R` and uses a manually created file defining which inversions span contigs in Hmel1 and Hmel2:
```
scripts/Figure3.R -i inversions/inversions.details.tsv -c inversions/contig_spanning_inversions.tsv -o Figure3
```


## Figures S6 and S7

Figures S6 and S7 were made from input file `inversions.gaps.tsv`, generated by `process_inversions.py`.
```r
library(ggplot2)
library(scales)
library(dplyr)
gaps<-read.delim("inversions/inversions.gaps.tsv")
maxchromlenmb<-19
pdf("FigureS6.pdf", width=16, height=6)
for (chrom in 1:21) {
  print(ggplot(gaps %>% filter(Chromosome==chrom), aes(GapLength, ymin=GapStart/1000000, ymax=GapEnd/1000000, colour=Species)) + 
          geom_linerange() + theme_bw() + coord_flip() + facet_grid(Species~.) +
          scale_y_continuous(limits=c(0,maxchromlenmb), breaks=seq(0,maxchromlenmb), minor_breaks=seq(0.5,maxchromlenmb)) + 
          scale_x_continuous(limits=c(0,1500000), labels=comma) +
          scale_colour_manual(values=c("#0072B2","#D55E00")) + 
          xlab("Gap Length (bp)") + ylab("Chromosome Position (Mb)") + ggtitle(paste("Chromosome",chrom)))
}
dev.off()

pdf("FigureS7.pdf",width=16, height=9)
ggplot(gaps,aes(GapLength)) + geom_histogram(binwidth=1000) + facet_grid(Species~.) + theme_bw() +
       scale_x_continuous(labels=function(x) x/1000,limits=c(1000,1500000), breaks=seq(0,1500000,100000)) + 
       xlab("Gap Length (kb)") + 
       scale_y_continuous(trans="log1p", breaks=c(0,1,2,5,10,20,50,100,200))
dev.off()
```


## Figure S8

Figure S8 was created using script `probability_of_missing_inversions.py`:

```
scripts/probability_of_missing_inversions.py -g inversions/inversions.gaps.tsv -o inversions.missing_probs.tsv
```

```r
library(ggplot2)
library(scales)
prob<-read.delim("inversions.missing_probs.tsv")
ggplot(prob, aes(Size, Probability, colour=Species)) + geom_line() + theme_bw() + scale_x_continuous(labels=comma, breaks=seq(0, 1500000, 100000)) + xlab("Inversion Size") + scale_y_continuous(limits=c(0,1), breaks=seq(0, 1, 0.1)) + ylab("Probability of detecting inversion") + scale_colour_manual(values=c("#0072B2", "#D55E00"))
ggsave("FigureS8.pdf")
```

## Figures S9 and S10

Read length histograms were generated from raw PacBio read FASTA files:
```
scripts/calculate_fasta_length_histogram.py -g 'cydno_females/*fasta'     -n cydno_females     > cydno_females.pacbio.length.hist
scripts/calculate_fasta_length_histogram.py -g 'cydno_males/*fasta'       -n cydno_males       > cydno_males.pacbio.length.hist
scripts/calculate_fasta_length_histogram.py -g 'melpomene_females/*fasta' -n melpomene_females > melpomene_females.pacbio.length.hist
scripts/calculate_fasta_length_histogram.py -g 'melpomene_males/*fasta'   -n melpomene_males   > melpomene_males.pacbio.length.hist
cat *.pacbio.length.hist > pacbio.length.histogram.tsv
```

Figure S9 was generated from these read length histograms:
```
library(ggplot2)
library(scales)
lengths<-read.delim("pacbio.length.histogram.tsv")
lengths$Sample<-revalue(lengths$Sample, c("cydno_females "="Cydno females", "cydno_males "="Cydno males", "melpomene_females "="Melpomene females", "melpomene_males "="Melpomene males"))
lengths$Sample<-factor(lengths$Sample, levels=c("Melpomene females","Cydno females","Melpomene males", "Cydno males"))
pdf("FigureS9.pdf",width=10,height=4)
ggplot(lengths, aes(Length,Reads,colour=Sample))+geom_line()+scale_x_sqrt(breaks=c(100,500,1000,2000,5000,6000,7500,10000,20000,30000,40000), labels=c(0.1,0.5,1,2,5,6,7.5,10,20,30,40))+theme_bw()+xlab("Read Length (kb)")+scale_colour_manual(values=c("#D55E00", "#0072B2", "#E69F00", "#56B4E9"))
dev.off()
```

Per-base coverage histograms were generated from primary alignments of PacBio reads to `Hmel2_ordered.fa` used to generate PBHoney candidate inversions:
```
samtools depth cydno_females.Hmel2_ordered.bwa.primary.bam     | scripts/calculate_pacbio_perbase_coverage.py > cydno_females.perbase.hist
samtools depth cydno_males.Hmel2_ordered.bwa.primary.bam       | scripts/calculate_pacbio_perbase_coverage.py > cydno_males.perbase.hist
samtools depth melpomene_females.Hmel2_ordered.bwa.primary.bam | scripts/calculate_pacbio_perbase_coverage.py > melpomene_females.perbase.hist
samtools depth melpomene_males.Hmel2_ordered.bwa.primary.bam   | scripts/calculate_pacbio_perbase_coverage.py > melpomene_males.perbase.hist
```

Figure S10 was generated from these per base coverage histograms:
```r
library(ggplot2)
library(scales)
library(dplyr)
cf<-read.delim("cydno_females.perbase.hist",header=FALSE)
mf<-read.delim("melpomene_females.perbase.hist",header=FALSE)
cm<-read.delim("cydno_males.perbase.hist",header=FALSE)
mm<-read.delim("melpomene_males.perbase.hist",header=FALSE)
cf$Sample<-"Cydno females"
cm$Sample<-"Cydno males"
mf$Sample<-"Melpomene females"
mm$Sample<-"Melpomene males"
depthhist<-bind_rows(cf, cm, mf, mm)
names(depthhist)<-c("Depth","Bases","Sample")
colours<-c("#0072B2", "#56B4E9", "#D55E00", "#E69F00")
pdf("FigureS10.pdf",width=10,height=6)
ggplot(depthhist, aes(Depth, Bases, fill=Sample, colour=Sample)) + geom_bar(stat="identity",position="dodge") + theme_bw() + 
       scale_x_continuous(limits=c(0,100),breaks=seq(0,100,5)) + scale_y_continuous(labels=comma,limits=c(0,22000000)) +
       scale_colour_manual(values=colours) + scale_fill_manual(values=colours)
dev.off()
```

## Population genetics statistics and figures S11-S17

Windows around inversions for calculation of population genetics statistics were calculated with script `make_inversion_windows.py`. This script creates a file for every inversion in a new output folder, which were merged together into one file.

```
scripts/make_inversion_windows.py -i inversions.details.tsv  -t genome/Hmel2.transitions.tsv -o inversion.windows
cat inversion.windows/*tsv > popgen/inversion.windows.tsv
```

Samples from Martin et al. 2013 were aligned to Hmel2 and VCFs parsed with [genomics_general](http://github.com/simonhmartin/genomics_general) (filtered `geno.csv` files provided in popgen folder):

```
for i in {0..21}; do (python ~/bin/genomics_general/VCF_processing/parseVCF.py -i martin2013.chr$i.vcf.gz --skipIndel --minQual 30 --gtf flag=DP min=10 max=200 --gtf flag=GQ min=20 siteTypes=SNP gtTypes=Het,HomAlt | python ~/bin/genomics_general/filterGenotypes.py -p ros ros.CJ2071,ros.CJ531,ros.CJ533,ros.CJ546 -p chi chi.CJ553,chi.CJ560,chi.CJ564,chi.CJ565 -p melG melG.CJ13435,melG.CJ9315,melG.CJ9316,melG.CJ9317 -p par par.JM371 -p ser ser.JM202 --minPopCalls 1,1,1,1,1 -o popgen/martin2013.chr$i.geno.csv &> martin2013.chr$i.parse.log &); done
```

Population genetics statistics were calculated as follows:
```
for i in {1..21}; do (for j in inversion.windows/inversion.windows.chr$i.*.tsv; do python ~/bin/genomics_general/popgenWindows.py --writeFailedWindows --windType predefined --windCoords $j -g popgen/martin2013.chr$i.geno.csv -p ros ros.CJ2071,ros.CJ531,ros.CJ533,ros.CJ546 -p chi chi.CJ553,chi.CJ560,chi.CJ564,chi.CJ565 -p melG melG.CJ13435,melG.CJ9315,melG.CJ9316,melG.CJ9317 -p par par.JM371 -p ser ser.JM202 -o ${j%tsv}popgen.csv -f phased -T 2 &> ${j%tsv}popgen.log; done &); done &
for i in {1..21}; do (for j in inversion.windows/inversion.windows.chr$i.*.tsv; do python ~/bin/genomics_general/ABBABABAwindows.py --writeFailedWindows --windType predefined --windCoords $j -g popgen/martin2013.chr$i.geno.csv -f phased -p P1 melG.CJ13435,melG.CJ9315,melG.CJ9316,melG.CJ9317 -p P2  ros.CJ2071,ros.CJ531,ros.CJ533,ros.CJ546 -p P3 chi.CJ553,chi.CJ560,chi.CJ564,chi.CJ565 -p O par.JM371,ser.JM202 -o ${j%tsv}abbababa.csv -T 2 &> ${j%tsv}abbababa.log; done &); done &
```

Output files were then merged together:
```
head -n1 inversion.windows/inversion.windows.chr1.1.popgen.csv  > popgen/inversion.windows.popgen.csv
grep -h ^chr inversion.windows/*.popgen.csv >> popgen/inversion.windows.popgen.csv
head -n1 inversion.windows/inversion.windows.chr1.1.abbababa.csv  > popgen/inversion.windows.abbababa.csv
grep -h ^chr inversion.windows/*.abbababa.csv >> popgen/inversion.windows.abbababa.csv
```

Figures S11-17 were then generated with script `make_inversion_plots.R`:

```
scripts/make_inversion_plots.R -i inversions/inversions.details.tsv -w popgen/inversion.windows.tsv -p popgen/inversion.windows.popgen.csv -a popgen/inversion.windows.abbababa.csv -s linkage_maps/linkage_map.crosses.map.tsv -c linkage_maps/linkage_map.species.cm.tsv -1 genome/Hmel1_chromosome_positions.tsv -o inversion.plots
```

# _H. erato_ comparison

The _H. erato_ genome assembly version 1 was downloaded from [LepBase](http://ensembl.lepbase.org/Heliconius_erato_v1/Info/Index) and aligned to `Hmel2_ordered.fa` with LAST:

```
lastdb -vcR11 Heratodb erato/Heliconius_erato_v1_-_scaffolds.fa
parallel-fasta -j 40 "lastal -v Heratodb | last-split" < Hmel2/Hmel2_ordered.fa > erato/Hmel2_ordered.Herato.maf
```

Script `Hmel2_Herato_maf.py` writes out chromosome positions for MAF alignments to `Hmel2_ordered.Herato.maf.tsv`, _H. erato_ scaffold lengths to `Herato_scaffolds.tsv`, and also the start positions of _H. erato_ scaffolds on _H. erato_ chromosomes to `Herato.chromosome_starts.tsv`.

```
scripts/Hmel2_Herato_maf.py -e genome/Herato_scaffolds.tsv -t genome/Hmel2.transitions.tsv -m genome/Hmel2_ordered.Herato.maf -o genome/Hmel2_ordered.Herato.maf.tsv > genome/Herato_scaffolds.tsv
```

_H. erato_ linkage maps were compiled into one file with script `compile_Herato_maps.py`:

```
scripts/compile_Herato_maps.py > genome/Herato_markers.tsv
```

Figure S18 was then created with script `Hmel2_Herato_dotplot.R`, which creates a PDF for each chromosome; these were merged together for the final figure:

```
scripts/Hmel2_Herato_dotplot.R  -a erato/Hmel2_ordered.Herato.maf.tsv -m genome/Hmel2_ordered.chromosome_starts.tsv -e genome/Herato.chromosome_starts.tsv -n linkage_map/linkage_map.species.cm.tsv -f genome/Herato_markers.tsv
```

