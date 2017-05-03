#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=FALSE))
library(ggplot2)
library(scales)
options(warn=-1)

# From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbpal <- c("#0072B2", "#D55E00", "#56B4E9", "#009E73", "#000000", "#E69F00", "#F0E442", "#CC79A7")
names(cbpal)<-c("blue","red","lightblue","green","black","orange","yellow","magenta")

maxchromlenmb<-19

parser <- ArgumentParser()
parser$add_argument("-i", "--inputstub", type='character', required=TRUE)
args <- parser$parse_args()

plot_cm<-function(cm, filesuffix) {
    pdf(paste(args$inputstub, filesuffix, 'marey', 'pdf',sep='.'),width=7.5, height=8.75)
    print(ggplot(cm, aes(End/1000000, cM, colour=Species, shape=Cross)) + geom_line() +
                 theme_bw(base_size=10) + theme(strip.background=element_blank()) +
                 scale_x_continuous(limits=c(0,maxchromlenmb), breaks=seq(0,maxchromlenmb,2), minor_breaks=seq(1,maxchromlenmb)) + 
                 facet_wrap(~Chromosome, nrow=7) + 
                 xlab("Chromosome Position (Mb)") + 
                 scale_colour_manual(values=c("#0072B2","#009E73","#D55E00"), guide=FALSE))
    dev.off()
}

species_cm<-read.delim(paste(args$inputstub, '.species.cm.tsv', sep=''))
species_cm$Cross<-'Species' # Dummy variable to allow generalisation with cross plot, has to be a factor so must be a string, not a number
out<-plot_cm(species_cm, 'species')

crosses_cm<-read.delim(paste(args$inputstub, '.crosses.cm.tsv', sep=''))
crosses_cm<-filter(crosses_cm, !Cross %in% c("C94", "C1"))
out<-plot_cm(crosses_cm, 'crosses')
