#!/usr/bin/env Rscript

# Author: John Davey johnomics@gmail.com

# Initialise
suppressPackageStartupMessages(library(argparse))
library(grid)
library(dplyr, warn.conflicts=FALSE)
library(ggplot2)
library(scales)
library(stringr)
suppressPackageStartupMessages(library(sqldf))
library(tcltk) # required for sqldf

# From http://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r#12135122
decimals <- function(x, k) format(round(x, k), nsmall=k)

# From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbpal <- c("#0072B2", "#D55E00", "#56B4E9", "#009E73", "#000000", "#E69F00", "#F0E442", "#CC79A7")
names(cbpal)<-c("blue","red","lightblue","green","black","orange","yellow","magenta")

sample_colours<-c("blue","lightblue","red","orange")
names(sample_colours)<-c("Cydno females", "Cydno males", "Melpomene females", "Melpomene males")

map_colours<-c(cbpal["green"], cbpal["yellow"])

stat_colours<-c("magenta","green","yellow")
stat_labels<-c(expression(F[ST]), expression(f[d]), expression(d[XY]))
stat_label_pos<-c(0.75,0.25,0.5)
names(stat_colours)<-names(stat_labels)<-names(stat_label_pos)<-c("Fst_ros_chi","fd","dxy_ros_chi")

species_colours<-c("red", "green", "blue")
names(species_colours)<-c("Melpomene","Hybrid","Cydno")

statuses<-c("Cydno_Both", "Cydno_SplitOnly",
            "Melpomene_Both", "Melpomene_SplitOnly",
            "Misassembly_Both", "Misassembly_SplitOnly", "Misassembly_Trio"
)

status_figure<-c("11","12","13","14","15","16","17")

status_species<-c(expression(italic("H. cydno")), expression(italic("H. cydno")),
                  expression(italic("H. melpomene")), expression(italic("H. melpomene")),
                  "Both species", "Both species", "Both species")

status_class<-c("Split reads and trio assembly", "Split reads only",
                "Split reads and trio assembly", "Split reads only", 
                "Split reads and trio assembly", "Split reads only", "Split reads in one species, trio assembly in both")

names(status_figure)<-names(status_species)<-names(status_class)<-statuses

options(width=200)

draw.species.map<-function(species, inv_snps, inv_cm, scaleStart, scaleEnd) {

    inv_cm <- filter(inv_cm, Species==species) %>%
              mutate(winprevcM = prevcM + cMFraction * (winStart-Start)/(End-Start),
                     wincM     = cM     - cMFraction * (1-((winEnd-Start)/(End-Start)))
              )

    if (nrow(inv_cm)>0) {
        grid.polyline(
            unit(c(inv_cm$winStart, inv_cm$winEnd), "native"),
            unit(c(inv_cm$winprevcM, inv_cm$wincM),"native"),
            id=rep(1:nrow(inv_cm),2),
            gp=gpar(col=cbpal[species_colours[species]], lwd=ifelse(inv_cm$cMFraction==0, 3, 1))
        )
    }


    inv_snps <- filter(inv_snps, Species==species) %>% arrange(Chromosome, ChromPosition)

    if (nrow(inv_snps)>0) {
        inv_snps <- sqldf("select s.ChromPosition, c.winStart, c.winEnd, c.winprevcM, c.wincM from inv_snps s inner join inv_cm c on (s.ChromPosition >= c.Start and s.ChromPosition <= c.End)")
        inv_snps <- mutate(inv_snps, cM = winprevcM + (wincM - winprevcM) * (ChromPosition-winStart)/(winEnd-winStart))
        grid.points(
            unit(inv_snps$ChromPosition, "native"),
            unit(inv_snps$cM, "native"),
            pch=16, size=unit(0.1,"inches"),
            gp=gpar(col=cbpal[species_colours[species]])
        )
    }
}

draw.map<-function(chromosome, scaleStart, scaleEnd, snps, cm, xlimits) {
    inv_snps <- filter(snps, Chromosome==chromosome, ChromPosition >= scaleStart, ChromPosition <= scaleEnd)
    inv_cm   <- filter(cm, Chromosome==chromosome) %>% group_by(Species) %>%
                mutate(prevcM = lag(cM, default=0)) %>% ungroup() %>%
                filter(Start <= scaleEnd, End >= scaleStart) %>%
                mutate(winStart = pmax(Start, scaleStart), winEnd = pmin(End, scaleEnd))

    mincM<-min(inv_cm$prevcM)
    maxcM<-max(inv_cm$cM)
    maxcM<-ifelse(mincM==maxcM, mincM+1,maxcM)
    pushViewport(viewport(0, 0.2, 1, 0.6, xscale=xlimits, yscale=c(mincM, maxcM), just=c("left","bottom")))

    grid.lines(unit(c(scaleStart,scaleStart),"native"), unit(c(mincM,maxcM),"native"), gp=gpar(lwd=2,col="grey60"))   # y-axis
    grid.lines(unit(c(scaleStart, scaleEnd),"native"), unit(c(mincM,mincM),"native"), gp=gpar(lwd=0.5, col="grey90")) # low  x-axis
    grid.lines(unit(c(scaleStart, scaleEnd),"native"), unit(c(maxcM,maxcM),"native"), gp=gpar(lwd=0.5, col="grey90")) # high x-axis

    tickStart<-scaleStart-(scaleEnd-scaleStart)*0.01
    grid.lines(unit(c(tickStart, scaleStart),"native"), unit(c(mincM,mincM),"native"), gp=gpar(lwd=2,col="grey60"))   # low  tick mark
    grid.lines(unit(c(tickStart, scaleStart),"native"), unit(c(maxcM,maxcM),"native"), gp=gpar(lwd=2,col="grey60"))   # high tick mark
    
    grid.text(paste(decimals(c(mincM, maxcM),3), "cM  "), unit(c(tickStart,tickStart),"native"), unit(c(mincM,maxcM),"native"), just=c("right","centre"),gp=gpar(fontsize=6))

    for (species in c("Cydno","Hybrid","Melpomene")) {
        draw.species.map(species, inv_snps, inv_cm, scaleStart, scaleEnd)
    }
    popViewport()
}

draw.key<-function() {
    grid.text(expression(italic("H. cydno")), 0.92, 0.97, just=c("centre","centre"), gp=gpar(fontsize=11,col=cbpal[sample_colours["Cydno females"]]))
    grid.text("females", 0.91, 0.95, just=c("right","centre"), gp=gpar(fontsize=8,col=cbpal[sample_colours["Cydno females"]]))
    grid.text("males", 0.93, 0.95, just=c("left","centre"), gp=gpar(fontsize=8,col=cbpal[sample_colours["Cydno males"]]))
    grid.text(expression(italic("H. melpomene")), 0.92, 0.90, just=c("centre","centre"), gp=gpar(fontsize=11,col=cbpal[sample_colours["Melpomene females"]]))
    grid.text("females", 0.91, 0.88, just=c("right","centre"), gp=gpar(fontsize=8,col=cbpal[sample_colours["Melpomene females"]]))
    grid.text("males", 0.93, 0.88, just=c("left","centre"), gp=gpar(fontsize=8,col=cbpal[sample_colours["Melpomene males"]]))
}

draw.hit.scaffold.candidates<-function(hits, species, ypos, arrows) {
    # Receives lines from one scaffold only (because draw.species.candidates groups by Label)
    scaffold<-unique(hits$Label)
    sex<-unique(hits$Sex)
    left<-min(hits$Start)
    right<-max(hits$End)
    if (length(sex) != 1) {
        print(hits)
        stop("More than one sex found!")
    }
    scf_ypos <- ypos + unique(hits$ScfPos)*0.4

    sample_colour <- cbpal[sample_colours[paste(species, sex)]]
    grid.text(paste(" ", scaffold), unit(max(hits$End),"native"), scf_ypos, just="left",gp=gpar(fontsize=5, col=sample_colour))

    pb_arrow<-NULL
    if (arrows) {
        pb_arrow<-arrow(ends=ifelse(hits$Dir==-1,"first","last"), angle=45, length=unit(0.03, "inches"))
    }
    grid.rect(unit(left,"native"),scf_ypos-0.03,unit(right-left,"native"),0.06, just=c("left","bottom"), gp=gpar(fill=sample_colour,lty="blank",alpha=0.2))
    grid.polyline(
        unit(c(hits$Start, hits$End), "native"),
        rep(scf_ypos+rep(c(0.01,-0.01),length.out=nrow(hits)), 2),
        id=rep(1:nrow(hits),2),
        gp=gpar(col=rep(sample_colour, 2), lwd=1),
        arrow=pb_arrow
    )
    hits
}

draw.pbhoney.candidates<-function(sample, ypos) {
    inv_ypos<-rep(rev(ypos+(1:nrow(sample))/nrow(sample)*0.4), 2)
    grid.polyline(
        unit(c(sample$Start, sample$End),"native"),
        inv_ypos,
        id=rep(1:nrow(sample), 2),
        gp=gpar(col=rep(sample$Colour, 2), lwd=3)
    )
    grid.rect(
        unit(sample$Start,"native"),
        -4.085,
        unit(sample$End-sample$Start,"native"),
        4.085+inv_ypos,
        just=c("left","bottom"),
        gp=gpar(fill=sample$Colour, lty="blank", alpha=0.02)
    )
    
    grid.text(paste(" ", sample$Reads), unit(sample$End,"native"),inv_ypos, just=c("left","centre"), gp=gpar(col=sample$Colour, fontsize=5))
}

draw.overlaps<-function(overlaps, species, ypos) {
    overlaps <- overlaps %>% filter(Species==species)
    if (nrow(overlaps) == 0) return()
    
    ov_ypos<-rep(rev(ypos+(1:nrow(overlaps))/nrow(overlaps)*0.4), 2)
    
    grid.polyline(
        unit(c(overlaps$ReadStart, overlaps$ReadEnd), "native"),
        ov_ypos,
        id=rep(1:nrow(overlaps), 2),
        gp=gpar(col=rep(cbpal[sample_colours[paste(overlaps$Species, overlaps$Sex)]], 2), lwd=0.5)
    )
}

draw.species.candidates<-function(lines, species, ypos, is.hit=FALSE, arrows=TRUE, overlaps=NULL) {
    if (!is.null(overlaps)) {
        draw.overlaps(overlaps, species, ypos)
    }

    sample<-filter(lines, Species==species)
    if (nrow(sample) == 0) return()

    sample$Colour<-cbpal[sample_colours[paste(species, sample$Sex)]]

    if (is.hit) {
        scaffolds<-length(unique(as.character(sample$Label)))
        sample$ScfPos<-as.integer(factor(sample$Label))/(scaffolds+1)
        sample %>% group_by(Label) %>% do(draw.hit.scaffold.candidates(., species, ypos, arrows))
    } else { # PBHoney
        draw.pbhoney.candidates(sample, ypos)
    }
}


draw.hittype<-function(inv, xlimits, hittype, hitlabel, ypos, rect_colour, is.hit=FALSE, arrows=TRUE, overlaps=NULL) {
    lines<-filter(inv, HitType==hittype)
    pushViewport(viewport(0,ypos, width=1, height=0.15, xscale=xlimits, just=c("left","bottom")))
    grid.rect(-0.1, 0, 0.175, 1, just=c("left","bottom"), gp=gpar(fill=rect_colour,lty="blank"))
    grid.text(hitlabel, 0, 0.5, rot=90, just=c("centre","top"), gp=gpar(fontsize=10))
    draw.species.candidates(lines, "Cydno", 0.5, is.hit, arrows, overlaps)
    draw.species.candidates(lines, "Melpomene", 0.0, is.hit, arrows, overlaps)
    popViewport()
}


draw.annotations<-function(inv) {
    feature_heights<-c(0.34,0.35,0.36,0.37,0.38,0.39)
    names(feature_heights)<-c("gene","rRNA","tRNA","mRNA","exon","CDS")
    grid.text(names(feature_heights), 0.925, feature_heights, just=c("left","centre"), gp=gpar(fontsize=5))
    gff<-filter(inv, HitType=="Feature")
    if (nrow(gff) == 0) return()
    grid.polyline(
        unit(c(gff$Start, gff$End), "native"),
        rep(feature_heights[as.character(gff$Label)],2),
        id=rep(1:nrow(gff),2),
        gp=gpar(col="brown", lwd=3, lineend="butt")
    )
}

draw.repeats<-function(inv) {
    rep<-filter(inv, HitType=="Repeat")
    if (nrow(rep) == 0) return()
    grid.polyline(
        unit(c(rep$Start, rep$End), "native"),
        rep(0.27, nrow(rep)*2),
        id=rep(1:nrow(rep), 2),
        gp=gpar(col="black", lwd=1, lineend="butt")
    )
}

draw.hmel2.contigs<-function(inv, scaleStart) {
    hmel2<-filter(inv, HitType=="Hmel2Part", Label != "Gap")
    grid.text("Genome",0,0.30,rot=90, just=c("centre","top"))
    grid.text("Hmel2  ", unit(scaleStart, "native"), 0.27, just=c("right","centre"), gp=gpar(fontsize=7))
    grid.polyline(
        unit(c(hmel2$Start, hmel2$End), "native"),
        rep(0.27, nrow(hmel2)*2),
        id=rep(1:nrow(hmel2),2),
        gp=gpar(col=rep(rep(c("snow3","snow4"),length.out=nrow(hmel2)),2), lwd=10, lineend="butt")
    )
    grid.text(hmel2$Label, unit(hmel2$Start+(hmel2$End-hmel2$Start)/2, "native"), rep(c(0.255,0.285),length.out=nrow(hmel2)), just=c("centre","centre"), gp=gpar(fontsize=5))
}

draw.hmel1.contigs<-function(hmel1, chromosome, scaleStart, scaleEnd) {
    hmel1 <- filter(hmel1, Chromosome==chromosome, ChromStart < scaleEnd, ChromEnd > scaleStart) %>% mutate(Start=pmax(ChromStart, scaleStart), End=pmin(ChromEnd, scaleEnd))
    grid.text("Hmel1  ", unit(scaleStart, "native"), 0.31, just=c("right","centre"), gp=gpar(fontsize=7))
    grid.polyline(
        unit(c(hmel1$Start, hmel1$End), "native"),
        rep(0.31, nrow(hmel1)*2),
        id=rep(1:nrow(hmel1),2),
        gp=gpar(col=rep(rep(c("snow3","snow4"),length.out=nrow(hmel1)),2), lwd=10, lineend="butt")
    )
    grid.text(hmel1$Hmel1Scaffold, unit(hmel1$Start+(hmel1$End-hmel1$Start)/2, "native"), rep(c(0.295,0.325),length.out=nrow(hmel1)), just=c("centre","centre"), gp=gpar(fontsize=5))
}

draw.physical.scale<-function(inv, scaleStart, scaleEnd) {
    scaleLen <- scaleEnd - scaleStart + 1
    xlimits<-c(scaleStart - scaleLen*0.1, scaleEnd + scaleLen*0.1)

    pushViewport(viewport(0.02, 0.02, width=0.96, height=0.96, xscale=xlimits, just=c("left","bottom")))

    pushViewport(viewport(0,0,1,0.8, xscale=xlimits, just=c("left","bottom"))) # 4/5ths scale to leave room for linkage map
    chromosome <- unique(inv$Chromosome)
    grid.text(paste("Chromosome", chromosome, "Position (bp)"),unit(scaleStart+scaleLen/2,"native"),0.19,just="centre",gp=gpar(fontsize=10))

    grid.lines(unit(c(scaleStart, scaleEnd), "native"), c(0.23,0.23), gp=gpar(lwd=3,col="grey70"))
    grid.text(format(round(c(scaleStart, scaleEnd)), big.mark=" "), unit(c(scaleStart, scaleEnd), "native"), 0.20, just="centre",gp=gpar(fontsize=8))
    grid.lines(unit(c(scaleStart,scaleStart), "native"), c(0.215,0.23), gp=gpar(lwd=3,col="grey70"))
    grid.lines(unit(c(scaleEnd,scaleEnd), "native"), c(0.215,0.23), gp=gpar(lwd=3,col="grey70"))

    group<-filter(inv, HitType=="Group")
    groupStart<-group$Start[1]
    groupEnd<-group$End[1]
    groupLen<-groupEnd - groupStart + 1
    grid.text(format(round(c(groupStart, groupEnd)), big.mark=" "), unit(c(groupStart, groupEnd), "native"), 0.2, just="centre",gp=gpar(fontsize=8))
    grid.lines(unit(c(groupStart,groupStart), "native"), c(0.215,0.23), gp=gpar(lwd=3,col="grey70"))
    grid.lines(unit(c(groupEnd,groupEnd), "native"), c(0.215,0.23), gp=gpar(lwd=3,col="grey70"))

    grid.lines(unit(c(groupStart, groupEnd), "native"), c(0.23, 0.23), gp=gpar(lwd=10, col="grey10"))
    grid.text(paste(format(round(groupLen), big.mark=" "), "bp"), unit((groupStart+groupEnd)/2, "native"), 0.23, just=c("centre","centre"), gp=gpar(fontsize=5, col="grey90"))

    xlimits
}

draw.stat<-function(winstats, stat, axis="left", statscale=1) {
    n<-nrow(winstats)
    winstats[[stat]]<-winstats[[stat]]/statscale
    invalid<-is.nan(winstats[[stat]]) | winstats[[stat]] < 0
    winstats$Shape<-ifelse(invalid,1,16)
    winstats[[stat]]<-ifelse(invalid,0,winstats[[stat]])
    grid.points(unit(winstats$midpoint, "native"),
                unit(winstats[[stat]], "native"),
                pch=winstats$Shape,
                gp=gpar(col=cbpal[stat_colours[stat]], cex=0.4)
    )
    winstats <- winstats %>% mutate(nextwin=lead(midpoint))
    winstats$nextstat <- c(winstats[[stat]][2:n], NA)
    grid.polyline(
        unit(c(winstats$midpoint, winstats$nextwin), "native"),
        unit(c(winstats[[stat]], winstats$nextstat), "native"),
        id=rep(1:n, 2),
        gp=gpar(col=cbpal[stat_colours[stat]])
    )
    grid.text(stat_labels[stat], ifelse(axis=="left",0,1.01), unit(stat_label_pos[stat],"native"),
              just=c(axis,"centre"), gp=gpar(col=cbpal[stat_colours[stat]], fontsize=6))
}

draw.popgen<-function(winstats) {
    scaleStart<-min(winstats$Start)
    scaleEnd<-max(winstats$End)
    scaleLen<-scaleEnd - scaleStart + 1

    invStart<-filter(winstats, window==0)$Start[1]
    invEnd<-filter(winstats, window==0)$End[1]

    pushViewport(viewport(0,0, width=1, height=0.17, xscale=c(scaleStart - scaleLen*0.05, scaleEnd + scaleLen*0.05), just=c("left","bottom")))

    # Plot x axis
    grid.lines(unit(c(scaleStart,scaleEnd), "native"), c(0.3,0.3), gp=gpar(lwd=2, col="grey60"))
    grid.text(format(c(scaleStart, invStart, invEnd, scaleEnd), big.mark=" "), unit(c(scaleStart, invStart, invEnd, scaleEnd), "native"), 0.2, just="centre",gp=gpar(fontsize=6))
    grid.polyline(
        unit(c(winstats$End, winstats$End),"native"),
        c(rep(0.25,nrow(winstats)), rep(0.3,nrow(winstats))),
        id=rep(1:nrow(winstats),2),
        gp=gpar(lwd=1,col="grey60")
    )
    grid.lines(unit(c(scaleStart,scaleStart), "native"), c(0.25,0.3), gp=gpar(lwd=2,col="grey60"))
    grid.lines(unit(c(scaleEnd,scaleEnd), "native"), c(0.25,0.3), gp=gpar(lwd=2,col="grey60"))

    # Plot inversion window axis
    grid.lines(unit(c(invStart, invEnd),"native"), c(0.3,0.3), gp=gpar(lwd=3, col="grey10"))
    
    # Plot stat
    pushViewport(viewport(0,0.3, width=1, height=0.7, xscale=c(scaleStart - scaleLen*0.05, scaleEnd + scaleLen*0.05), yscale=c(0,1), just=c("left","bottom")))
    grid.lines(unit(c(scaleStart,scaleEnd), "native"), unit(c(0.5,0.5), "native"), gp=gpar(lwd=0.5, col="grey90"))
    grid.lines(unit(c(scaleStart,scaleEnd), "native"), unit(c(1,1), "native"), gp=gpar(lwd=0.5, col="grey90"))
    winstats$midpoint<-winstats$Start + (winstats$End-winstats$Start)/2
    draw.stat(winstats, "dxy_ros_chi", "right", 0.1)
    draw.stat(winstats, "Fst_ros_chi", "left")
    draw.stat(winstats, "fd", "left")
    popViewport()

    # Plot left y axis
    grid.lines(unit(c(scaleStart, scaleStart), "native"), c(0.3,1), gp=gpar(lwd=2,col="grey60"))
    tickstart<-scaleStart-scaleLen*0.01
    grid.lines(unit(c(tickstart, scaleStart),"native"), c(0.3,   0.3  ), gp=gpar(lwd=2,col="grey60"))
    grid.lines(unit(c(tickstart, scaleStart),"native"), c(0.65,  0.65 ), gp=gpar(lwd=1,col="grey60"))
    grid.lines(unit(c(tickstart, scaleStart),"native"), c(1, 1),           gp=gpar(lwd=2,col="grey60"))
    grid.text(c("0 ", "0.5 ", "1 "), unit(rep(tickstart, 5), "native"), c(0.3, 0.65, 1), just=c("right","centre"), gp=gpar(fontsize=4))

    # Plot right y axis
    grid.lines(unit(c(scaleEnd, scaleEnd), "native"), c(0.3,1), gp=gpar(lwd=2,col="grey60"))
    tickend<-scaleEnd+scaleLen*0.01
    grid.lines(unit(c(scaleEnd, tickend),"native"), c(0.3,   0.3  ), gp=gpar(lwd=2,col="grey60"))
    grid.lines(unit(c(scaleEnd, tickend),"native"), c(0.65,  0.65 ), gp=gpar(lwd=1,col="grey60"))
    grid.lines(unit(c(scaleEnd, tickend),"native"), c(1, 1),           gp=gpar(lwd=2,col="grey60"))
    grid.text(c(" 0", " 0.05", " 0.1"), unit(rep(tickend, 5), "native"), c(0.3, 0.65, 1), just=c("left","centre"), gp=gpar(fontsize=4))

    
    # Plot site counts
    winstats$siteratio<-winstats$sites/max(winstats$sites)
    winstats$sitewidth<-(invEnd-invStart)*winstats$siteratio
    grid.polyline(
        unit(c(winstats$midpoint - winstats$sitewidth/2, winstats$midpoint + winstats$sitewidth/2),"native"),
        c(rep(0.02,nrow(winstats)*2)),
        id=rep(1:nrow(winstats),2),
        gp=gpar(lwd=2, col="grey70")
    )
    grid.text(format(winstats$sites, big.mark=' '), unit(winstats$midpoint,"native"), rep(0.05,nrow(winstats)), just=c("centre","bottom"), gp=gpar(fontsize=5, col="grey70"))
    grid.text("Sites",unit(scaleEnd + scaleLen*0.02, "native"), 0.04, just=c("centre","centre"),gp=gpar(fontsize=6,col="grey70"))

    popViewport()

}


plot_page<-function(inv, winstats, snps, cm, overlaps, hmel1) {
    grid.newpage()

    group<-filter(inv, HitType=="Group")
    groupStart<-group$Start[1]
    groupEnd<-group$End[1]
    groupLen<-groupEnd-groupStart+1
    scaleStart<-groupStart-groupLen/2
    scaleEnd<-groupEnd+groupLen/2
    chromosome <- unique(inv$Chromosome)
    inv %>% filter(Start < scaleEnd, End > scaleStart) %>% mutate(Start=pmax(Start,scaleStart),End=pmin(End,scaleEnd)) -> inv

    grid.text(paste("Figure S", status_figure[group$GroupStatus], '.', group$order, sep=''), 0.02, 0.98, just=c("left","top"), gp=gpar(fontsize=10))
    grid.text(status_species[group$GroupStatus], 0.5, 0.98, just=c("centre","top"), gp=gpar(fontsize=10))
    grid.text(status_class[group$GroupStatus], 0.98, 0.98, just=c("right","top"), gp=gpar(fontsize=10))
    xlimits<-draw.physical.scale(inv, scaleStart, scaleEnd)
    
    if (!is.null(overlaps)) {
        overlaps <- overlaps %>% filter(Chromosome==chromosome, ReadStart < scaleEnd & ReadEnd > scaleStart)
    }
    draw.hittype(inv, xlimits, "Inverted", "Trio\nInverted", 0.7, "grey88", is.hit=TRUE)
    draw.hittype(inv, xlimits, "Spanning", "Trio\nSpanning", 0.55, "grey92",  is.hit=TRUE)
    draw.hittype(inv, xlimits, "Edge",    "Trio\nEdges",    0.4, "grey96", is.hit=TRUE)
    draw.hittype(inv, xlimits, "PBHoney",  "PBHoney\nCandidates",       0.85, "grey84", FALSE, overlaps=overlaps)

    draw.annotations(inv)
    
    draw.hmel2.contigs(inv, scaleStart)

    draw.repeats(inv)

    draw.hmel1.contigs(hmel1, chromosome, scaleStart, scaleEnd)
    
    draw.key()

    draw.popgen(winstats %>% filter(candidate==as.character(group$candidate[1])))
    
    popViewport() # 4/5ths scale
    pushViewport(viewport(0,0.8,1,0.2,xscale=xlimits, just=c("left","bottom"))) # 1/5th scale for map row
    draw.map(chromosome, scaleStart, scaleEnd, snps, cm, xlimits)
    popViewport() # Map row
    popViewport() # Page

    inv
}

plot_group<-function(group, outputstub, winstats, snps, cm, overlaps, hmel1) {

  groupstatus<-group$GroupStatus[1]
  message(paste("Plotting",groupstatus))
  pdf(paste(outputstub, groupstatus, "pdf", sep='.'), width=7, height=9)
  
  grouplen<-filter(group, HitType=="Group") %>% select(GroupStatus, GroupID, GroupLength=Length) %>%
            arrange(desc(GroupLength), GroupID) %>% mutate(order = 1:length(GroupID))

  group<-merge(group, grouplen) %>% group_by(GroupID) %>% mutate(Label=factor(Label, rev(levels(Label))))

  for (i in 1:max(group$order)) {
      plot_page(filter(group, order==i), winstats, snps, cm, overlaps, hmel1)
  }
  dev.off()
  
  data.frame()
}

load.popgen<-function(winfile, popgenfile, abbababafile) {
    windows<-read.delim(winfile, sep=' ', header=FALSE)
    names(windows)<-c("scaffold","start","end","candidate","window","iteration")
    popgen<-read.csv(popgenfile)
    abbababa<-read.csv(abbababafile)

    winstats<-full_join(full_join(popgen,windows),abbababa)
    winstats$Chromosome<-as.numeric(str_extract(winstats$scaffold, "[0-9]+"))
    winstats<-rename(winstats, Start=start, End=end)
    
    winstats
}

args<-commandArgs(trailingOnly=T)

parser <- ArgumentParser()
parser$add_argument('-i', '--inversions', type='character', required=TRUE)
parser$add_argument('-o', '--outputstub', type='character', default="test")
parser$add_argument('-w', '--windows', type='character', required=TRUE)
parser$add_argument('-p', '--popgen', type='character', required=TRUE)
parser$add_argument('-a', '--abbababa', type='character', required=TRUE)
parser$add_argument('-s', '--snps', type='character', required=TRUE)
parser$add_argument('-c', '--cm', type='character', required=TRUE)
parser$add_argument('-r', '--read_overlaps', type='character')
parser$add_argument('-1', '--hmel1', type='character', required=TRUE)
args <- parser$parse_args()


message("Loading inversions")
groups<-read.delim(args$inversions)
groups$HitType<-factor(groups$HitType, c("Group","PBHoney","Tentative","Inverted","Spanning","Edge", "Feature", "Hmel2Part", "Repeat","Reject","Simulated"))

mergecats<-function(x) {
    x<-as.character(x)
    x<-ifelse(x=="Cydno_PBHoney" | x=="Cydno_Single", "Cydno_SplitOnly", x)
    x<-ifelse(x=="Melpomene_PBHoney" | x=="Melpomene_Single", "Melpomene_SplitOnly", x)
    x<-ifelse(x=="Misassembly_PBHoney", "Misassembly_SplitOnly", x)
    factor(x)
}

groups$Status<-mergecats(groups$Status)
groups$GroupStatus<-mergecats(groups$GroupStatus)

groups$Status<-as.character(groups$Status)
levels(groups$Label)<-gsub(" ", "\n", levels(groups$Label))

message("Loading popgen")
winstats<-load.popgen(args$windows, args$popgen, args$abbababa)
groups<-full_join(groups, winstats %>% filter(window==0, iteration==0) %>% select(Chromosome, Start, End, candidate))

message("Loading maps")
snps<- read.delim(args$snps) %>% select(Species, Chromosome, ChromPosition) %>% unique()
cm<-read.delim(args$cm)

overlaps<-NULL
if (!is.null(args$read_overlaps)) {
    message("Loading read overlaps")
    overlaps<-read.delim(args$read_overlaps)
}

message("Loading Hmel1 positions")
hmel1<-read.delim(args$hmel1, header=FALSE)
names(hmel1)<-c("Chromosome", "ChromStart", "ChromEnd", "Hmel2Scaffold", "Hmel2Start", "Hmel2End", "Hmel1Scaffold", "Hmel1Start", "Hmel1End")

groups %>% filter(GroupStatus != "Reject") %>% group_by(GroupStatus) %>% do(plot_group(., args$outputstub, winstats, snps, cm, overlaps, hmel1)) -> nullout
