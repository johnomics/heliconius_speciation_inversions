#!/usr/bin/env Rscript

# Author: John Davey johnomics@gmail.com

# Initialise
library(grid)
library(dplyr, warn.conflicts=FALSE)
suppressPackageStartupMessages(library(argparse))
library(parallel)
suppressPackageStartupMessages(library(sqldf))

options(width=200)

args<-commandArgs(trailingOnly=T)

genetic.colours<-c(rgb(244,165,130,max=255),rgb(5,113,176,max=255))

species<-c(expression(italic("H. cydno")), expression(italic("H. cydno")), expression(italic("H. cydno")), expression(italic("H. melpomene")), expression(italic("H. melpomene")),  expression(italic("H. melpomene")), "Hybrids")
crossnames<-c("Cross 1", "Cross 2", "Cross 3", "Cross 1", "Cross 2", "Cross 3", "")
names(crossnames)<-names(species)<-c("CYDA","CYDB","CYDE","MEL1","MEL3","MEL4","Hybrid")


draw.physical.scale<-function(maxchrsize) {
    pushViewport(viewport(0.05,0,width=0.95,height=0.05, xscale=c(0, maxchrsize), just=c("left","bottom"))) # Physical scale
    mb.onemil.tick<-seq(0,maxchrsize,1000000)
    chrmb<-ceiling(maxchrsize/1000000)
    grid.text("Chromosome length (Megabases)",unit(maxchrsize/2,"native"),0.1,just="centre",gp=gpar(fontsize=10))
    grid.text(sprintf("%2d", 0:chrmb), unit(mb.onemil.tick, "native"), 0.6, just="centre",gp=gpar(fontsize=8))
    grid.lines(unit(c(0,maxchrsize), "native"), c(1,1), gp=gpar(lwd=3, col="grey",lineend="round"))
    grid.polyline(
        unit(c(mb.onemil.tick,mb.onemil.tick),"native"),
        c(rep(0.8,chrmb),rep(1,chrmb)),
        id=rep(1:chrmb,2),
        gp=gpar(col="grey",lwd=2,lineend="round")
    )
    popViewport()
}

draw.physical.map<-function(scaffolds, maternals, maxchrsize, chr, cross) {
    
    scfcol<-apply(scaffolds, 1, function(x) {if (x["PoolType"]=="ok") "#009E73" else if (x["PoolType"]=="orient") "#E69F00" else if (as.numeric(x["PoolID"]) %% 2 == 0) "indianred1" else "indianred"})

    grid.polyline(
        unit(c(scaffolds$ChromStart,scaffolds$ChromEnd),"native"),
        unit(c(rep(0.8,nrow(scaffolds)),rep(0.8,nrow(scaffolds))),"native"),
        id=rep(1:nrow(scaffolds),2),
        gp=gpar(col=scfcol,lwd=13,lineend="butt")
    )
    
    grid.text(paste(cross, scaffolds$Hmel2Refined),unit(scaffolds$ChromStart,"native"),0.8,just=c("centre","top"),rot=90,gp=gpar(fontsize=1))
    
    if (nrow(maternals) == 0) return()
    grid.polyline(
        unit(c(maternals$ChromPosition,maternals$ChromPosition),"native"),
        rep(c(0.87,0.91),each=nrow(maternals)),
        id=rep(1:nrow(maternals),2),
        gp=gpar(lwd=0.1, lineend="butt", col="grey25")
    )
}

draw.linkage.map<-function(markers, ypos, maxcmlength) {
    geneticvp<-viewport(xscale=c(0, maxcmlength))
    if (nrow(markers)==0) return(geneticvp)

    pushViewport(geneticvp)


    maxcm <- max(markers$cM)
    
    grid.lines(unit(c(0,maxcm),"native"),c(0.2,0.2),gp=gpar(col=rgb(141,160,203,max=255),lwd=3,lineend="round"))
    grid.text(sprintf("%.3f",markers$cM),unit(markers$cM,"native"),0.14,just=c("right","centre"),rot=90,gp=gpar(fontsize=3))
    grid.polyline(
        unit(c(markers$cM,markers$cM),"native"),
        c(rep(0.15,nrow(markers)),rep(0.25,nrow(markers))),
        id=rep(1:nrow(markers),2),
        gp=gpar(col=markers$Colour,lwd=1,lineend="round")
    )

    popViewport()

    return(geneticvp)
}


plot.snps<-function(markers, snps, geneticvp) {
    if (nrow(snps)==0) return()
    snps<-merge(snps, data.frame(cM=markers$cM,Colour=markers$Colour), by="cM")
    snps$Colour<-as.character(snps$Colour)

    grid.polyline(
        unit(c(snps$ChromPosition,snps$ChromPosition),"native"),
        c(rep(0.62,nrow(snps)),rep(0.72,nrow(snps))),
        id=rep(1:nrow(snps),2),
        gp=gpar(col=snps$Colour,lwd=0.5,lineend="round")
    )

    apply (snps,1,
        function(x) {
            names(x)<-names(snps)
            pushViewport(geneticvp)
            grid.move.to(unit(x["cM"],"native"),0.25)
            popViewport()
            grid.line.to(unit(x["ChromPosition"],"native"),0.62,gp=gpar(col=x["Colour"],lwd=0.25,lty="dashed"))
        }
    )
}

plotcross<-function(chr, i, cross, proportion, markers, snps, maternals, scaffolds, maxchrsize, maxcmlength) {
    cross<-strsplit(cross,"\\.")[[1]][1]
    message(paste("Plotting chromosome", chr, "cross", cross))
    ypos<-1-i*proportion

    pushViewport(viewport(0,ypos,width=1,height=proportion,just=c("left","bottom")))
    
    if (cross %in% names(species)) {
        grid.text(species[cross], 0.045, 0.8, just=c("right", "centre"), gp=gpar(fontsize=6))
        grid.text(crossnames[cross], 0.045, 0.6, just=c("right","centre"), gp=gpar(fontsize=6))
    }
    if (cross != "Hybrid") {
        grid.text(cross, 0.045, 0.4, just=c("right","centre"), gp=gpar(fontsize=6))
    }

    pushViewport(viewport(0.05,0,width=0.95,height=1,xscale=c(0,maxchrsize),just=c("left","bottom")))
    draw.physical.map(scaffolds, maternals, maxchrsize, chr, cross)

    markers$Colour<-rep_len(genetic.colours,nrow(markers))
    geneticvp<-draw.linkage.map(markers, ypos, maxcmlength)

    plot.snps(markers, snps, geneticvp)
    
    popViewport()
    popViewport()

}

plotpage<-function(chr, snps, maternals, transitions, maxchrsize, maxcmlength, output) {
    outputstub<-ifelse(is.null(output),"test",output)
    outfilename<-paste(outputstub, ".", chr, ".pdf", sep="")

    crosses <- sort(as.character(unique(snps$Cross)))
    proportion <- 1/(length(crosses))

    pdf(outfilename, width=11.69, height=max(8.27, length(crosses)*1.1))

    grid.newpage()

    pushViewport(viewport(0.01, 0.01, width=0.96, height=0.96, just=c("left","bottom"))) #Page

    draw.physical.scale(maxchrsize)
    grid.text(paste("Chromosome",chr), 0.5, 1, just=c("centre","top"))

    pushViewport(viewport(0,0.06,width=1,height=0.92,just=c("left","bottom"))) #Plotting region

    for (i in 1:length(crosses)) {
        cross <- crosses[i]
        markers<-distinct(select(snps[snps$Chromosome==chr & snps$Cross==cross,], Chromosome, cM))
        markers<-markers[with(markers, order(Chromosome, cM)),]
        plotcross(chr, i, cross, proportion, markers,
                  snps[snps$Chromosome==chr & snps$Cross==cross,], 
                  maternals[maternals$Chromosome==chr & maternals$Cross==cross,],
                  transitions[transitions$Chromosome==chr,], maxchrsize, maxcmlength)
    }
    popViewport() #Plotting region
    popViewport() #Page
    dev.off()
}


load.data<-function(filename) {
    read.delim(filename, colClasses=c(Pattern="character"))
}

parser <- ArgumentParser()
parser$add_argument("-t", "--transitions", type='character', required=TRUE)
parser$add_argument("-o", "--outputstub", type='character')
parser$add_argument("-s", "--snps", type='character', required=TRUE)
parser$add_argument("-m", "--maternals", type='character', required=TRUE)
parser$add_argument('-a', "--hybridsnps", type='character', required=TRUE)
parser$add_argument('-l', "--linkagegroup", type='integer') # Single chromosome for output
args <- parser$parse_args()

message("Loading map")
read.delim(args$transitions)->transitions
transitions <- transitions %>% filter(Chromosome > 0) %>% group_by(Chromosome) %>% mutate(ChromLength=max(ChromEnd))

message("Loading SNPs")
data<-mclapply(c(args$snps, args$maternals, args$hybridsnps), load.data, mc.cores=4)
snps<-data[[1]]
maternals<-data[[2]]
hybridsnps<-data[[3]]

maternals<-subset(maternals, MarkerType == 2)

snps<-filter(snps, Species != "Hybrid")
snps<-rbind(snps, hybridsnps)

maxchrsize  <- ceiling(max(transitions$ChromLength)/1000000)*1000000+1
maxcmlength <- max(snps$cM)

crosses<-unique(snps$Cross)

chrs<-sort(unique(transitions$Chromosome))
if(!is.null(args$linkagegroup)) chrs<-args$linkagegroup

nullout<-mcmapply(plotpage, chrs, MoreArgs=list(snps, maternals, transitions, maxchrsize, maxcmlength, args$outputstub), mc.cores=4)

if (!is.null(warnings())) {
    warnings()
}
quit()
