#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
library(parallel)
library(dplyr, warn.conflicts=FALSE)
library(tidyr)
library(ggplot2)
library(scales)
suppressPackageStartupMessages(library(doMC))
library(foreach)

# From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbpal <- c("#0072B2", "#D55E00", "#56B4E9", "#009E73", "#000000", "#E69F00", "#F0E442", "#CC79A7")
names(cbpal)<-c("blue","red","lightblue","green","black","orange","yellow","magenta")

cols<-c("#0072B2","#009E73","#D55E00")

offspring<-c(297,335,331,95,77,125,111,122,102,170,88,68,5)
names(offspring)<-c("Cydno","Melpomene","Hybrid","CYDA","CYDB","CYDE","MEL1","MEL3","MEL4","C6","C26","C29","C1")
maxchromlenmb<-19

chromlengths<-c(17206585, 9045316, 10541528, 9662098, 9908586, 14054175, 14308859, 9320449, 8708747, 17965481, 11759272, 16327298, 18127314, 9174305, 10235750, 10083215, 14773299, 16803890, 16399344, 14871695, 13359691)

get_num_offspring<-function(species) {
    species<-unique(species)
    num_offspring<-offspring[species]
    names(num_offspring)<-NULL
    return(num_offspring)
}

recfrac<-function(RecStart, RecEnd, WinStart, WinEnd, resample=FALSE) {
  
  if (resample) {
    recsample <- sample(length(RecStart), replace=TRUE)
    RecStart <- RecStart[recsample]
    RecEnd   <- RecEnd[recsample]
  }
  
  RecLength = RecEnd - RecStart + 1

  FracStart <- pmax(RecStart, WinStart)
  FracEnd   <- pmin(RecEnd,   WinEnd  )
  
  recs <- ifelse(FracStart < FracEnd, (FracEnd-FracStart+1)/RecLength, 0)
  
  return(round(sum(recs),1))
}

get_rec_window <- function(recombinations, WinStart, WinEnd, windowsize) {
    recombinations %>% group_by(Group, Chromosome) %>% filter(Start < End, WinStart <= End, Start < WinEnd) %>%
    summarise(windowpos=WinStart + windowsize/2,
              recs=recfrac(Start, End, WinStart, WinEnd),
              cm=round(recs/get_num_offspring(Group) * 100 * windowsize/1000000, 3))
}

calculate_rec_rates <- function(recombinations) {
    recwindows <- foreach(WinStart=seq(1,19000000,args$step), .combine='bind_rows', .multicombine=TRUE, .maxcombine=1000) %do% {
                        WinEnd <- WinStart + args$windowsize - 1
                        get_rec_window(recombinations, WinStart, WinEnd, args$windowsize)
                  }
    lastwin <- foreach(chr=1:21, .combine='bind_rows') %do% {
        get_rec_window(recombinations %>% filter(Chromosome==chr), chromlengths[chr] - args$windowsize + 1, chromlengths[chr], args$windowsize)
    }
    lastwin <- ungroup(lastwin) %>% complete(Group, nesting(Chromosome, windowpos), fill=list(recs=0, cm=0))

    recwindows <- recwindows %>% complete(Group, Chromosome, windowpos, fill=list(recs=0, cm=0)) %>%
                                 bind_rows(lastwin) %>% 
                                 filter(windowpos <= chromlengths[Chromosome]-args$windowsize/2 + 1)
    recwindows
}

plot_rec_rate<-function(recwindows, outputtype, outputstub) {
    pdf(paste(outputstub, ".recombination_rate.", outputtype, ".pdf", sep=""), width=16, height=6)
    for (chrom in 1:21) {
        recplot <- ggplot(recwindows %>% filter(Chromosome==chrom), aes(windowpos/1000000, cm)) + geom_line(aes(colour=Species, shape=Cross))
        maxy <- max(recwindows$cm)
        if ('cm.qu025' %in% names(recwindows)) {
            recplot <- recplot + geom_ribbon(aes(ymin=cm.qu025, ymax=cm.qu975, fill=Species), alpha=0.1, linetype="blank") + scale_fill_manual(values=cols)
            maxy <- max(c(recwindows$cm, recwindows$cm.qu975))
        }
        
        print(recplot + scale_x_continuous(limits=c(0,maxchromlenmb), breaks=seq(0,maxchromlenmb), minor_breaks=seq(0.5,maxchromlenmb)) + 
                        scale_y_continuous(limits=c(0,maxy)) + theme_bw() + scale_colour_manual(values=cols, guide=FALSE) +
                        ylab("Recombination rate (cM/Mb)") + xlab("Chromosome Position (Mb)") + ggtitle(paste("Chromosome",chrom)))
    }
    dev.off()
    return()
}

parser <- ArgumentParser()
parser$add_argument("-r", "--recombinations", type='character', required=TRUE)
parser$add_argument("-t", "--threads", type='integer', default=1)
parser$add_argument("-i", "--iterations", type='integer', default=1)
parser$add_argument("-o", "--outputstub", type='character', default="test")
parser$add_argument("-w", "--windowsize", type='integer', default=1000000)
parser$add_argument("-s", "--step", type='integer', default=100000)
args <- parser$parse_args()

start.time <- Sys.time()
registerDoMC(args$threads)

message("Loading recombinations")
recombinations<-read.delim(args$recombinations,stringsAsFactors=FALSE)

message("Calculating cross recombination rates")
crosses.recwindows <- calculate_rec_rates(recombinations %>% rename(Group=Cross)) %>% rename(Cross=Group) %>%
                      inner_join(recombinations %>% select(Species, Cross) %>% distinct(), by="Cross") %>%
                      arrange(Species, Cross, Chromosome, windowpos)
write.table(crosses.recwindows, file=paste(args$outputstub, ".windows.crosses.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

message("Plotting cross recombination rates")
plot_rec_rate(crosses.recwindows, 'crosses', args$outputstub)


message("Calculating species recombination rates")
species.recwindows <- calculate_rec_rates(recombinations %>% rename(Group=Species)) %>% rename(Species=Group)

if (args$iterations > 1) {
    message("Calculating species bootstrap confidence intervals")
    species.recs <- recombinations %>% rename(Group=Species) %>% select(-Cross)
    species.offspring<- species.recs %>% select (Group, Offspring) %>% group_by(Group) %>% unique()
    bootstrap.recwindows <- foreach (i=1:args$iterations, .combine='bind_rows', .inorder=FALSE, .multicombine=TRUE, .maxcombine=1000) %dopar% {
        if (i %% max(round(args$iterations/10),1) == 0) message(paste("Calculating iteration",i))
        offspring.i <- species.offspring %>% group_by(Group) %>% mutate(Offspring=sample(Offspring, replace=TRUE))
        recs.i <- inner_join(offspring.i, species.recs, by=c("Group","Offspring"))
        recwindows.i <- calculate_rec_rates(recs.i)
        recwindows.i$Iteration <- i
        recwindows.i
    }
    species.ranges <- bootstrap.recwindows %>% rename(Species=Group) %>% group_by(Species, Chromosome, windowpos) %>%
                      summarise(recs.qu025=quantile(recs, 0.025), recs.qu975=quantile(recs, 0.975), cm.qu025=quantile(cm, 0.025), cm.qu975=quantile(cm, 0.975))
    species.recwindows<-merge(species.recwindows, species.ranges) %>% arrange(Species, Chromosome, windowpos)
}

write.table(species.recwindows, file=paste(args$outputstub, ".windows.species.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

message("Plotting species recombination rates")
species.recwindows$Cross<-species.recwindows$Species
plot_rec_rate(species.recwindows, 'species', args$outputstub)

end.time<-Sys.time()
run.time<-end.time - start.time
message(paste("Run time:", as.integer(run.time), attributes(run.time)$units))