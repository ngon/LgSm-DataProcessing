######## Auxillary code to generate plots ###############

## Working directory on Natalia's computer PalmerLab-Atlas
setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/figures/gwas")


## First, read in chrLengths from the mm10 reference.
## In CRI: 'group/palmer-lab/reference_genomes/mouse/chrLengths_mm10_chrnames.txt'
chrLens    <- read.table("./dataFiles/chrLengths_mm10.txt")$V2
chrLens    <- c(0, cumsum(chrLens*1.0))

labPos     <- read.table("./dataFiles/chrLengths_mm10.txt")$V2/2.0
labPos     <- chrLens[1:19]+labPos


traits <- c("cM.act1.t", "cM.startle", "cM.sc8.t", "cM.sc1.t")
,"wild.binary", "act5.t")


titles <- c(#"Change in preference for meth paired side",
            #"Time spent on meth paired side on day 8 (5 min)",
            #"Time spent on meth paired side on day 8 (30 min)",
            "Locomotor response to novelty (distance traveled on Day 1)",
            #"Locomotor response to meth (Day 2 activity)",
#             "Day 4 activity (1 mg/kg meth)",
#             "Day 5 activity (saline)",
#             "Day 8 activity (saline)",
            #"Locomotor sensitization to meth",
            "Acoustic startle reflex (120 dB)",
            "Exploratory activity Day 8 (number of side changes)",
            "Exploratory activity Day 1 (number of side changes)",
            "Wildness (number of escapes during CPP test)",
            "Locomotor activity (distance traveled on Day 5)")
#             "Prepulse inhibition (73 dB)",
#             "Prepulse inhibition (76 dB)",
#             "Prepulse inhibition (82 dB)",
            #"Habituation to acoustic startle",
            #"Blood glucose levels after a 4h fast",
#            "Tail length (cm)"
            )

outfiles <- paste0(traits, ".manhattan.png")


for (index in 1:4) {
        print(paste("Starting trait", traits[index]))
        plotManhattan(traits[index], outfile=outfiles[index], title=titles[index])
        print(paste("Done with trait", traits[index]))
}

plotManhattan(trait="wild.binary", outfile="wildbinary.manhattan.png", signifLine=7,
              title="Wildness (number of escapes during CPP test)")

plotManhattan(trait="act5.t", outfile="act5t.manhattan.png", signifLine=6.25,
              title="Locomotor activity (distance traveled on Day 5)")
# colors
# light grey = #aab0be
# coral = #ff6040
# steelblue = #63B8FF

plotManhattan <- function(trait, outfile=NULL,
                          oddcolor="#63B8FF", evencolor="#aab0be",
                          signifLine=6, signifColor="#000000", title=NULL) {
        if (is.null(outfile)) {
                outfile <- paste0(trait, ".manhattan.png")
        }
        if (is.null(title)) {
                title = trait
        }
        png(file=outfile, height=600, width=1400, bg="transparent")
        pvals <- c()
        positions <- c()
        cols <- c()

        for (chrom in 1:19) {
                filename <- paste0(trait, ".chr", chrom, ".assoc.txt")
                print(paste("Starting chromosome", chrom))
                output <- read.table(filename, header=T, as.is=TRUE)
                output$chr <- as.numeric(sub("chr", "", output$chr))
                pvals <- c(pvals, -log10(output$p_lrt))
                positions <- c(positions, chrLens[output$chr]+output$ps)
                if (chrom %% 2) {
                        cols <- c(cols, rep(oddcolor, nrow(output)))
                } else {
                        cols <- c(cols, rep(evencolor, nrow(output)))
                }
                print(paste("Done with chromosome", chrom))
        }
        plot(positions, pvals, col=cols, pch=19, cex=0.78, main=title, ylab='-log10(p-value)', xlab='Position', axes=F, frame.plot=T, cex.lab=1.7, cex.main=1.8, font.main=1, ylim =c(0,7) )
        axis(2, cex.axis=1.5)
        axis(1, at=labPos, labels=1:19, tick=F, cex.axis=1.5)
        abline(h=signifLine, col=signifColor, lty=2)
        dev.off()
}

__________________________________________

### zoom plot for cpp.diff - ch4

# cppdiff.ch4 <- read.table("cpp.diff.chr4.assoc.txt", sep="\t", header=T )
# test<-cppdiff.ch4[order(cppdiff.ch4$p_lrt),]
# plotZoom(trait="cpp.diff", chrom=4, start=89803883, end=114867413)
#
# ppi6.ch8 <- read.table("ppi6.chr8.assoc.txt", sep="\t", header=T )
# test<-ppi6.ch8[order(ppi6.ch8$ps),]
# plotZoom(trait="ppi6", chrom=8, start=69822096, end=95001306)
#
# ppi12.ch8 <- read.table("ppi12.chr8.assoc.txt", sep="\t", header=T )
# test<-ppi12.ch8[order(ppi12.ch8$ps),]
# plotZoom(trait="ppi12", chrom=8, start=69822096, end=95001306)
#
# startle17 <- read.table("startle.chr17.assoc.txt", sep="\t", header=T )
# test<-startle17[order(startle17$ps),]
# plotZoom(trait="startle", chrom=17, start=18998527, end=35964216)
#
# startle7 <- read.table("startle.chr7.assoc.txt", sep="\t", header=T )
# test<-startle7[order(startle7$ps),]
# plotZoom(trait="startle", chrom=7, start=75001287, end=84997578)
#
# cpp5min.7 <- read.table("cpp8.1.chr7.assoc.txt", sep="\t", header=T )
# test<-cpp5min.7[order(cpp5min.7$ps),]
# plotZoom(trait="cpp8.1", chrom=7, start=75001287, end=99997517)
#
# act1t7 <- read.table("act1.t.chr7.assoc.txt", sep="\t", header=T )
# test<-act1t7[order(act1t7$ps),]
# plotZoom(trait="act1.t", chrom=7, start=75001287, end=99997517)
#
# act1t17 <- read.table("act1.t.chr17.assoc.txt", sep="\t", header=T )
# test<-act1t17[order(act1t17$ps),]
# plotZoom(trait="act1.t", chrom=17, start=, end=)
#
# act2t7 <- read.table("act2.t.chr7.assoc.txt", sep="\t", header=T )
# test<-act2t7[order(act2t7$ps),]
# plotZoom(trait="act2.t", chrom=7, start=30000032, end=45847045)
#
# act2t17 <- read.table("act2.t.chr17.assoc.txt", sep="\t", header=T )
# test<-act2t17[order(act2t17$ps),]
# plotZoom(trait="act2.t", chrom=17, start=75999608, end= 85035207)

act2t4 <- read.table("act2.t.chr4.assoc.txt", sep="\t", header=T )
test<-act2t4[order(act2t4$ps),]
plotZoom(trait="act2.t", chrom=4, start=127158422, end=134999602)

act5t4 <- read.table("act5.t.chr4.assoc.txt", sep="\t", header=T )
test<-act5t4[order(act5t4$ps),]
plotZoom(trait="act5.t", chrom=4, start=127158422, end=134999602)


hab <- read.table("habituation.chr4.assoc.txt", sep="\t", header=T )
test<-hab[order(hab$ps),]
head(test)
tail(test)

plotZoom(trait="habituation", chrom=4, start=6063986, end=156255616)


plotZoom <- function(trait, chrom, start, end, outfile=NULL, dotColor="#aab0be", signifLine=5, signifColor="#000000", title=NULL, scaling=1e6, scalename="(Mb)") {
    if (is.null(outfile)) {
        outfile <- paste0(trait, ".chr", chrom, ".", as.character(start), ".", as.character(end), ".zoom.jpg")
    }
    if (is.null(title)) {
        title = trait
    }
    jpeg(file=outfile, quality=95)
    filename <- paste0(trait, ".chr", chrom, ".assoc.txt")
    print(paste("Starting chromosome", chrom))
    output <- read.table(filename, header=T)
    pvals <- -log10(output$p_lrt)
    positions <- output$ps/scaling
    print(paste("Done with chromosome", chrom))
    ### Code to plot
    start = start/scaling
    end = end/scaling
    plot(positions, pvals, col=dotColor, pch=19, cex=0.75, main=title, ylab='-log10(p-value)', xlab=paste('Position', scalename), axes=F, frame.plot=T, cex.lab=1.1, xlim=c(start, end))
    axis(2, cex.axis=1.2)
    axis(1, cex.axis=1.2)
    abline(h=signifLine, col=signifColor, lty=2)
    dev.off()
}


3158013
94657175


