#!/usr/bin/env Rscript

# Author: John Davey johnomics@gmail.com

# Initialise
suppressPackageStartupMessages(library(argparse))
library(plyr)
library(dplyr, warn.conflicts=FALSE)
library(ggplot2)
library(scales)



# From http://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r#12135122
decimals <- function(x, k) format(round(x, k), nsmall=k)

# From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbpal <- c("#0072B2", "#D55E00", "#56B4E9", "#009E73", "#000000", "#E69F00", "#F0E442", "#CC79A7","gray70", "gray40")
names(cbpal)<-c("blue","red","lightblue","green","black","orange","yellow","magenta","gray70", "gray40")

status_figure<-c("11","12","13","14","15","16","17")
status_colours<-c("blue","lightblue","red","orange","gray40","gray70","gray70")
status_cbcol<-cbpal[status_colours]
names(status_figure)<-names(status_cbcol)<-c("Cydno_Both", "Cydno_SplitOnly", 
                        "Melpomene_Both", "Melpomene_SplitOnly",
                        "Misassembly_Both", "Misassembly_SplitOnly", "Misassembly_Trio")

# http://stackoverflow.com/questions/33524669/labeling-outliers-of-boxplots-in-r
is_outlier <- function(x) {
  return(x > quantile(x, 0.75) + 1.5 * IQR(x))
}

options(width=300)

args<-commandArgs(trailingOnly=T)

parser <- ArgumentParser()
parser$add_argument('-i', '--inversions', type='character', required=TRUE)
parser$add_argument('-o', '--output', type='character', default="test")
parser$add_argument('-c', '--contigspanning', type='character', required=TRUE)
args <- parser$parse_args()

details<-read.delim(args$inversions)

mergecats<-function(x) {
    x<-as.character(x)
    x<-ifelse(x=="Cydno_PBHoney" | x=="Cydno_Single", "Cydno_SplitOnly", x)
    x<-ifelse(x=="Melpomene_PBHoney" | x=="Melpomene_Single", "Melpomene_SplitOnly", x)
    x<-ifelse(x=="Misassembly_PBHoney", "Misassembly_SplitOnly", x)
    factor(x)
}

details$Status<-mergecats(details$Status)
details$GroupStatus<-mergecats(details$GroupStatus)

details<-filter(details, GroupID > 0)

contigs<-read.delim(args$contigspanning)
details<-merge(details, contigs)

grouplen<-filter(details, HitType=="Group") %>% group_by(GroupStatus) %>%
          select(GroupStatus, GroupID, GroupLength=Length) %>%
          mutate(outlier = ifelse(is_outlier(GroupLength), GroupLength, as.numeric(NA))) %>%
          arrange(desc(GroupLength), GroupID) %>% mutate(order = 1:length(GroupID))

details<-merge(details, grouplen)

details$Class  <- sapply(as.character(details$GroupStatus), function(x) strsplit(x, '_')[[1]][1])
details$Status <- sapply(as.character(details$GroupStatus), function(x) strsplit(x, '_')[[1]][2])

ClassLabels<-c("H. cydno", "H. melpomene", "Both species")
names(ClassLabels)<-c("Cydno", "Melpomene", "Misassembly")
details$ClassLabels<-factor(ClassLabels[details$Class], levels=c("H. cydno", "H. melpomene", "Both species"))

StatusLabels<-c("Split reads and\ntrio assembly", "Split reads only", "Split reads in\none species only,\ntrio assembly in both")
names(StatusLabels)<-c("Both", "SplitOnly", "Trio")
details$StatusLabels<-StatusLabels[details$Status]

pdf(paste(args$output, "pdf", sep='.'), width=7.5, height=5.55)
plotdetails <- details %>% filter(HitType=="Group") %>% arrange(GroupStatus, order)
ggplot(plotdetails, aes(Status, Length, fill=GroupStatus))  + 
       geom_boxplot(outlier.shape=NA, varwidth=TRUE) +
       geom_point(aes(y=outlier, shape=Hmel1 | Hmel2), size=2) +
       geom_text(aes(label=paste('   S', status_figure[GroupStatus], '.', order, '   ', sep=''), y=outlier, hjust=rep(c("left","right"),length.out=length(order))), size=2) + 
       facet_wrap(~ClassLabels, scales="free_x") +
       theme_bw(base_size=9) + theme(legend.position="bottom") + 
       scale_x_discrete(labels=StatusLabels) + scale_fill_manual(values=status_cbcol, guide=FALSE) + 
       scale_y_log10(labels=comma, breaks=c(1000, 5000, 10000, 50000, 100000, 500000)) + scale_shape_manual(values=c(16, 4), guide=FALSE) +
       xlab("Inversion Evidence") + ylab("Inversion Length (bp)")

dev.off()