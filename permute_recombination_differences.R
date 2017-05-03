#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr, warn.conflicts=FALSE))
library(tidyr)
library(ggplot2)
library(scales)
suppressPackageStartupMessages(library(doMC))
library(foreach)

offspring<-c(297,335,331)
names(offspring)<-c("Cydno","Melpomene","Hybrid")
maxchromlenmb<-19

chromlengths<-c(17206585, 9045316, 10541528, 9662098, 9908586, 14054175, 14308859, 9320449, 8708747, 17965481, 11759272, 16327298, 18127314, 9174305, 10235750, 10083215, 14773299, 16803890, 16399344, 14871695, 13359691)

na.zero <- function (x) {
    x[is.na(x)] <- 0
    return(x)
}

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
    recombinations %>% group_by(Species, Chromosome) %>%
    filter(Start < End, WinStart <= End, Start < WinEnd) %>%
    summarise(WindowPos=WinStart + windowsize/2,
        recs=recfrac(Start, End, WinStart, WinEnd),
        cm=round(recs/get_num_offspring(Species)*100 * windowsize/1000000,3))
}


calculate_differences <- function(recombinations) {
    recwindows <- foreach(WinStart=seq(1,19000000,args$step), .combine='bind_rows', .multicombine=TRUE, .maxcombine=1000) %do% {
                        WinEnd <- WinStart + args$windowsize - 1
                        get_rec_window(recombinations, WinStart, WinEnd, args$windowsize)
                   }

    lastwin <- foreach (chr=1:21, .combine='bind_rows') %do% {
        get_rec_window(recombinations %>% filter(Chromosome==chr), chromlengths[chr] - args$windowsize + 1, chromlengths[chr], args$windowsize)
    }
    lastwin <- ungroup(lastwin) %>% complete(Species, nesting(Chromosome, WindowPos), fill=list(recs=0, cm=0))
    
    recwindows <- recwindows %>% complete(Species, Chromosome, WindowPos, fill=list(recs=0, cm=0)) %>%
                    bind_rows(lastwin) %>% 
                    filter(WindowPos <= chromlengths[Chromosome] - args$windowsize/2 + 1) %>%
                    select(Species, Chromosome, WindowPos, cm) %>%
                    spread(Species, cm) %>% mutate_each(funs(na.zero)) %>%
                    mutate(MelCyd = abs(Melpomene-Cydno), MelHyb=abs(Melpomene-Hybrid), CydHyb=abs(Cydno-Hybrid)) %>%
                    gather(Comparison, Difference, MelCyd, MelHyb, CydHyb)
    recwindows
}


parser <- ArgumentParser()
parser$add_argument("-r", "--recombinations", type='character', required=TRUE)
parser$add_argument("-t", "--threads", type='integer', default=1)
parser$add_argument("-p", "--permutations", type='integer', default=1)
parser$add_argument("-o", "--outputstub", type='character', default="test")
parser$add_argument("-w", "--windowsize", type='integer', default=1000000)
parser$add_argument("-s", "--step", type='integer', default=100000)
parser$add_argument("-d", "--dump", action="store_true", default=FALSE)
parser$add_argument("-l", "--load", type='character', default='')
args <- parser$parse_args()

start.time <- Sys.time()
registerDoMC(args$threads)

recombinations<-read.delim(args$recombinations, stringsAsFactors=FALSE)
if (is.null(recombinations$Iteration)) recombinations$Iteration<-0

offspring_names <- recombinations %>% group_by(Species) %>% select(Species, Offspring) %>% do(Offspring=unique(.$Offspring)) %>% unnest(Offspring)

message(paste("Calculating cM differences for", as.integer(args$windowsize), "bp windows in", as.integer(args$step), "bp steps using", args$threads, "threads"))

message("Calculating true differences")
truewindows<-calculate_differences(recombinations)
setnames(truewindows, old=c('Comparison', 'Difference'), new=c('TrueComparison', 'TrueDifference'))

if (args$load != '') {
	message(paste("Loading permutations from", args$load))
	permwindows<-read.delim(args$load)
} else {
	message("Generating permutations")
	permrecs <- recombinations
	permwindows<-foreach(perm = 1:args$permutations, .combine='bind_rows', .inorder=FALSE, .multicombine=TRUE, .maxcombine=10000) %dopar% {
	    if (perm %% max(round(args$permutations/10),1) == 0) message(paste("Calculating permutation",perm))
	    permrecs$Species <- NULL
	    permoffspring <- data.frame(Species=sample(offspring_names$Species), Offspring=offspring_names$Offspring)
	    permrecs <- merge(permrecs, permoffspring) %>% group_by(Species, Chromosome) %>% arrange(Start)
	    permrecs$Species <- as.character(permrecs$Species)
	    permwin<-calculate_differences(permrecs) %>% select(-Cydno, -Melpomene, -Hybrid)
	    permwin$Iteration <- perm
	    permwin
	}
}

if (args$dump) {
	write.table(permwindows, paste(args$outputstub, ".windowdump.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
} else {
	message("Calculating summary stats")
	pvals <- permwindows %>% group_by(Chromosome, WindowPos, Comparison) %>% left_join(truewindows %>% select(Chromosome, WindowPos, TrueComparison, TrueDifference), by=c("Chromosome","WindowPos","Comparison"="TrueComparison")) %>%
	         summarise(TrueDifference = round(unique(TrueDifference),3), p = (sum(TrueDifference<Difference)+1)/(n()+1),
	                   Qu95= round(quantile(Difference,0.95),2), Qu99= round(quantile(Difference,0.99),2), Max=round(max(Difference),3))
	
	write.table(pvals, file=paste(args$outputstub, ".cmdifferences.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)
}

truewindows <- truewindows %>% select(-TrueComparison, -TrueDifference) %>% distinct(Chromosome, WindowPos, Melpomene, Cydno, Hybrid) %>% gather(Species, cM, Melpomene, Cydno, Hybrid) %>% arrange(Chromosome, WindowPos, Species)
write.table(truewindows, file=paste(args$outputstub, ".cm.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

end.time<-Sys.time()
run.time<-end.time - start.time
message(paste("Run time:", as.integer(run.time), attributes(run.time)$units))
