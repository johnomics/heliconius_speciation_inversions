#!/usr/bin/env Rscript

# Author: John Davey johnomics@gmail.com

# Initialise
suppressPackageStartupMessages(library(argparse))
library(grid)
library(dplyr, warn.conflicts=FALSE)
library(scales)
suppressPackageStartupMessages(library(doMC))
library(foreach)

# From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbpal <- c("#0072B2", "#D55E00", "#56B4E9", "#009E73", "#000000", "#E69F00", "#F0E442", "#CC79A7")
names(cbpal)<-c("blue","red","lightblue","green","black","orange","yellow","magenta")

melchromlen<-c(17206585, 9045316, 10541528, 9662098, 9908586, 14054175, 14308859, 9320449, 8708747, 17965481, 11759272, 16327298, 18127314, 9174305, 10235750, 10083215, 14773299, 16803890, 16399344, 14871695, 13359691)

erachromlen<-c(22325789, 13638029, 16642927, 15722921, 13990926, 19430166, 19097038, 14025409, 13164215, 23964976, 17158922, 21359818, 23854410, 14030661, 17148182, 15262899, 20434887, 21837394, 22828410, 19449475, 17478929)

options(width=200)

args<-commandArgs(trailingOnly=T)

parser <- ArgumentParser()
parser$add_argument('-a', '--alignments', type='character', required=TRUE)
parser$add_argument('-m', '--melstarts', type='character', required=TRUE)
parser$add_argument('-e', '--eratostarts', type='character', required=TRUE)
parser$add_argument('-n', '--melmap', type='character', required=TRUE)
parser$add_argument('-f', '--eramap', type='character', required=TRUE)
parser$add_argument('-t', '--threads', type='integer', default=4)
args <- parser$parse_args()

message("Loading alignments")
alignments<-read.delim(args$alignments, stringsAsFactors=FALSE)

melstarts<-read.delim(args$melstarts)
erastarts<-read.delim(args$eratostarts)

melmap<-read.delim(args$melmap) %>% filter(Species=="Melpomene") %>% group_by(Chromosome, cM) %>% summarise(Start=min(Start), End=max(End))
eramap<-read.delim(args$eramap)

registerDoMC(args$threads)

null<-foreach (chrom=1:21)  %dopar% {
    melmapchrom<-filter(melmap, Chromosome==chrom)
    eramapchrom<-filter(eramap, Chromosome==chrom)
    melchromstarts<-filter(melstarts, Chromosome==chrom)
    erachromstarts<-filter(erastarts, Chromosome==chrom)
    chromalign <- filter(alignments, Hmel2_Chromosome == chrom, Herato_Chromosome == chrom)
    
    pdf(paste("Hmel2_Herato_dotplot", chrom, "pdf",sep='.'), width=ceiling(erachromlen[chrom]/1000000), height=ceiling(melchromlen[chrom]/1000000))
    message(paste("Plotting chromosome", chrom))
    grid.newpage()
    pushViewport(viewport(x=0.02, y=0.02, w=0.96, h=0.96, just=c("left","bottom")))
    
    grid.text(chrom, 0, 1, just=c("left","top"),gp=gpar(font=2, fontsize=24))
    # Melpomene scaffold names
    pushViewport(viewport(x=0, y=0.1, w=0.1, h=0.8, yscale=c(melchromlen[chrom],1), just=c("left","bottom")))
    grid.text(melchromstarts$Scaffold, x=0.9, y=unit(melchromstarts$ChromStart + melchromstarts$Length/2, "native"), just=c("right","centre"))
    popViewport()

    # Erato scaffold names
    pushViewport(viewport(x=0.1, y=0.9, w=0.8, h=0.1, xscale=c(1,erachromlen[chrom]), just=c("left","bottom")))
    grid.text(erachromstarts$Scaffold, unit(erachromstarts$ChromStart + erachromstarts$Length/2, "native"), 0.1, just=c("right","centre"), rot=270)
    popViewport()
    
    # Main window
    pushViewport(viewport(x=0.1, y=0.1, w=0.8, h=0.8, xscale=c(1,erachromlen[chrom]), yscale=c(melchromlen[chrom],1), just=c("left","bottom")))
    grid.polygon(unit(c(0, erachromlen[chrom], erachromlen[chrom], 0),"native"), unit(c(0, 0, melchromlen[chrom], melchromlen[chrom]), "native"))
    
    mel_cm_col<-rep(c("gray70", "white"), length.out=nrow(melmapchrom))
    grid.rect(
        rep(0, nrow(melmapchrom)),
        unit(c(melmapchrom$Start), "native"),
        rep(1, nrow(melmapchrom)),
        unit(c(melmapchrom$End-melmapchrom$Start+1), "native"),
        just=c("left","bottom"),
        gp=gpar(fill=mel_cm_col, lty="blank", alpha=0.1)
    )

    era_cm_col<-rep(c("gray70", "white"), length.out=nrow(eramapchrom))
    grid.rect(
        unit(c(eramapchrom$Start), "native"),
        rep(0, nrow(eramapchrom)),
        unit(c(eramapchrom$End-eramapchrom$Start+1), "native"),
        rep(1, nrow(eramapchrom)),
        just=c("left","bottom"),
        gp=gpar(fill=era_cm_col, lty="blank", alpha=0.1)
    )

    
    grid.polyline(
        rep(c(0,1), each=nrow(melchromstarts)),
        unit(rep(melchromstarts$ChromStart, 2), "native"),
        id=rep(1:nrow(melchromstarts), 2),
        gp=gpar(col="gray80")
    )
    
    grid.polyline(
        unit(rep(erachromstarts$ChromStart, 2), "native"),
        rep(c(0,1), each=nrow(erachromstarts)),
        id=rep(1:nrow(erachromstarts), 2),
        gp=gpar(col="gray80")
    )
        
    for (i in 1:nrow(chromalign)) {
        a<-chromalign[i,]

        melstart<-a$Hmel2_ChromStart
        melend  <-a$Hmel2_ChromEnd
        if (a$Hmel2_Orientation == '-') {
            melstart<-a$Hmel2_ChromEnd
            melend  <-a$Hmel2_ChromStart
        }

        erastart<-a$Herato_ChromStart
        eraend  <-a$Herato_ChromEnd
        if (a$Herato_Orientation == '-') {
            erastart<-a$Herato_ChromEnd
            eraend  <-a$Herato_ChromStart
        }

        acol<-ifelse(a$Hmel2_Orientation == a$Herato_Orientation, cbpal['red'], cbpal['blue'])
        
        grid.lines(unit(c(erastart, eraend), "native"), unit(c(melstart, melend), "native"), gp=gpar(col=acol))
    }
    
    grid.polyline(
        unit(c(eramapchrom$Start, eramapchrom$End), "native"),
        rep(-0.005, nrow(eramapchrom)*2),
        id=rep(1:nrow(eramapchrom), 2),
        gp=gpar(col=c("black", "gray50"), lineend="butt", lwd=10)
    )

    grid.polyline(
        rep(1.005, nrow(melmapchrom)*2),
        unit(c(melmapchrom$Start, melmapchrom$End), "native"),
        id=rep(1:nrow(melmapchrom), 2),
        gp=gpar(col=c("black", "gray50"), lineend="butt", lwd=10)
    )
    
    popViewport()
    
    # Erato axis
    pushViewport(viewport(x=0.1, y=0, w=0.8, h=0.1, xscale=c(1,erachromlen[chrom]), just=c("left","bottom")))
    grid.lines(unit(c(0, erachromlen[chrom]),"native"), c(0.7,0.7), gp=gpar(lwd=3, col="gray70"))
    erato_mb<-seq(0, erachromlen[chrom], 1000000)
    grid.polyline(
        unit(rep(erato_mb, 2), "native"),
        rep(c(0.5, 0.7), each=length(erato_mb)),
        id=rep(1:length(erato_mb), 2),
        gp=gpar(lwd=2, col="gray70")
    )
    grid.text(erato_mb/1000000, unit(erato_mb, "native"), rep(0.4, length(erato_mb)), gp=gpar(fontsize=10))
    grid.text("H. erato Chromosome Position (Mb)", 0.5, 0.15, gp=gpar(fontsize=12))
    popViewport()
    
    # Melpomene axis
    pushViewport(viewport(x=0.9, y=0.1, w=0.1, h=0.8, yscale=c(melchromlen[chrom],1), just=c("left","bottom")))
    grid.lines(c(0.3,0.3), unit(c(0, melchromlen[chrom]),"native"), gp=gpar(lwd=3, col="gray70"))
    mel_mb<-seq(0, melchromlen[chrom], 1000000)
    grid.polyline(
        rep(c(0.3, 0.5), each=length(mel_mb)),
        unit(rep(mel_mb, 2), "native"),
        id=rep(1:length(mel_mb), 2),
        gp=gpar(lwd=2, col="gray70")
    )
    grid.text(mel_mb/1000000, rep(0.6, length(mel_mb)), unit(mel_mb, "native"), rot=270, gp=gpar(fontsize=10))
    grid.text("H. melpomene Chromosome Position (Mb)", 0.85, 0.5, rot=270, gp=gpar(fontsize=12))
    popViewport()
    dev.off()
}
