#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
library(parallel)
library(dplyr, warn.conflicts=FALSE)
library(tidyr)
library(ggplot2)
library(grid)
library(scales)
suppressPackageStartupMessages(library(doMC))
library(foreach)

cols<-c("#0072B2","#009E73","#D55E00")

hybridfamilies<-c("CYDA","CYDB","CYDE","MEL1","MEL3","MEL4","C8","C18","C23","C10","C24","C14","C11","C20","C32","C37","C36","C51","C38","C72","C70","C71","C94","C3")
hybridparents<-c("CYDA","CYDB","CYDE","MEL1","MEL3","MEL4","C6","C6","C6","C6","C6","C6","C6","C6","C26","C26","C26","C26","C26","C29","C29","C29","C29","C1")
hybridcross<-data.frame(ParentCross=hybridparents,Cross=hybridfamilies, stringsAsFactors=FALSE)

maxchromlenmb<-19

chromlengths<-c(17206585, 9045316, 10541528, 9662098, 9908586, 14054175, 14308859, 9320449, 8708747, 17965481, 11759272, 16327298, 18127314, 9174305, 10235750, 10083215, 14773299, 16803890, 16399344, 14871695, 13359691)

get_snp_density <- function(snps, WinStart, WinEnd, windowsize) {
    snps %>% group_by(Group, Chromosome) %>% filter(WinStart <= ChromPosition, ChromPosition <= WinEnd) %>%
    summarise(windowpos=WinStart + windowsize/2,
              snps=length(unique(ChromPosition))
             )
}

calculate_snp_densities <- function(snps) {
    windows <- foreach(WinStart=seq(1,19000000,args$step), .combine='bind_rows', .multicombine=TRUE, .maxcombine=1000) %do% {
                        WinEnd <- WinStart + args$windowsize - 1
                        get_snp_density(snps, WinStart, WinEnd, args$windowsize)
                  }
    windows <- windows %>% complete(Group, Chromosome, windowpos, fill=list(snps=0)) %>%
                                 filter(windowpos <= chromlengths[Chromosome]-args$windowsize/2 + 1)
    windows
}

plot_snp_density<-function(windows, outputtype, outputstub) {
    pdf(paste(outputstub, ".snp_density.", outputtype, ".pdf", sep=""), width=16, height=12)
    xscale <- scale_x_continuous(limits=c(0,maxchromlenmb), breaks=seq(0,maxchromlenmb), minor_breaks=seq(0.5,maxchromlenmb))
    for (chrom in 1:21) {
        winplot<-ggplot(windows %>% filter(Chromosome==chrom) %>% rename(GCPC=GC) %>% mutate(PstI="               ", GC="               "), aes(x=windowpos/1000000)) +
                 xscale + theme_bw()
        snpplot <- winplot + geom_line(aes(y=snps, colour=Species, shape=ParentCross)) + scale_colour_manual(values=cols) +
                             ylab("SNPs per Megabase\n") + ggtitle(paste("Chromosome", chrom)) +
                             scale_y_continuous(limits=c(0, max(windows$snps))) + theme(axis.title.x=element_blank())
        resplot <- winplot + geom_line(aes(y=Ressites, colour=PstI, shape=ParentCross)) +
                             scale_colour_manual(values="black") +
                             ylab("PstI sites per Megabase\n") + theme(axis.title.x=element_blank()) +
                             scale_y_continuous(limits=c(0, max(windows$Ressites)))
        gcplot <- winplot + geom_line(aes(y=GCPC, colour=GC, shape=ParentCross)) +
                             scale_colour_manual(values="black") +
                             ylab("GC content per Megabase (%)\n") + xlab("Chromosome Position (Mb)") +
                             scale_y_continuous(limits=c(30, 40))


        
        grid.newpage()
        grid.draw(rbind(ggplotGrob(snpplot), ggplotGrob(resplot), ggplotGrob(gcplot), size="last"))
    }
    dev.off()
    return()
}

options(width=200)
parser <- ArgumentParser()
parser$add_argument("-r", "--snps", type='character', required=TRUE)
parser$add_argument("-t", "--threads", type='integer', default=1)
parser$add_argument("-o", "--outputstub", type='character', default="test")
parser$add_argument("-w", "--windowsize", type='integer', default=1000000)
parser$add_argument("-s", "--step", type='integer', default=100000)
parser$add_argument("-f", "--gcressite", type='character', required=TRUE)
parser$add_argument("-c", "--cm", type="character", required=TRUE)

args <- parser$parse_args()

start.time <- Sys.time()
registerDoMC(args$threads)

message("Loading SNPs")
snps<-read.delim(args$snps, stringsAsFactors=FALSE)
snps<-inner_join(snps, hybridcross, by="Cross")

message("Loading GC and ressite stats")
gc<-read.delim(args$gcressite, stringsAsFactors=FALSE)

message("Loading recombination rates")
cm<-read.delim(args$cm, stringsAsFactors=FALSE) %>% select(Species, Chromosome, windowpos, cm)

message("Calculating species SNP densities")
species.windows <- calculate_snp_densities(snps %>% rename(Group=Species)) %>% rename(Species=Group) %>%
                      inner_join(snps %>% select(Species, ParentCross) %>% distinct(), by="Species") %>%
                      arrange(Species, Chromosome, windowpos)
species.windows<-inner_join(species.windows, gc, by=c("Chromosome","windowpos"="Position"))

plot_snp_density(species.windows, 'species', args$outputstub)

species.cor<-inner_join(species.windows, cm, by=c("Species","Chromosome","windowpos")) %>% select(-ParentCross) %>% distinct(Species, Chromosome, windowpos, snps, GC, Ressites, cm, .keep_all=TRUE)
print(species.cor %>% group_by(Species) %>% summarise(
    PstICor=cor.test(snps, Ressites, alternative="greater")$estimate,
    PstIp  =cor.test(snps, Ressites, alternative="greater")$p.value,
    PstIsig=cor.test(snps, Ressites, alternative="greater")$p.value<0.05/3,
    
    cMCor=cor.test(snps, cm, alternative="greater")$estimate,
    cMp  =cor.test(snps, cm, alternative="greater")$p.value,
    cMsig=cor.test(snps, cm, alternative="greater")$p.value<0.05/3,
    
    GCCor=cor.test(GC, Ressites, alternative="greater")$estimate,
    GCp  =cor.test(GC, Ressites, alternative="greater")$p.value,
    GCsig=cor.test(GC, Ressites, alternative="greater")$p.value<0.05/3) %>% as.data.frame())

write.table(species.cor, file=paste(args$outputstub, ".snp_density.windows.species.tsv", sep=""), sep="\t", row.names=FALSE, 
quote=FALSE)


message("Calculating cross SNP densities")
crosses.windows <- calculate_snp_densities(snps %>% rename(Group=ParentCross)) %>% rename(ParentCross=Group) %>%
                      inner_join(snps %>% select(Species, ParentCross) %>% distinct(), by="ParentCross") %>%
                      arrange(Species, ParentCross, Chromosome, windowpos)
crosses.windows<-inner_join(crosses.windows, gc, by=c("Chromosome"="Chromosome","windowpos"="Position"))
plot_snp_density(crosses.windows, 'crosses', args$outputstub)
write.table(crosses.windows, file=paste(args$outputstub, ".snp_density.windows.crosses.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

