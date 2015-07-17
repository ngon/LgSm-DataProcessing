#######  PLOTTING COMMANDS  --------------------------------------------------------------
library("ggplot2")
source("./scripts/multiplot.R")
load("./dataFiles/snpDensity.allData.RData")
load("./dataFiles/ibd.intervals.RData")
#source("/group/palmer-lab/AIL/LgSm-DataProcessing/scripts/multiplot.R")
#load("/group/palmer-lab/AIL/LgSm-DataProcessing/dataFiles/snpDensity.allData.RData")
#load("/group/palmer-lab/AIL/LgSm-DataProcessing/dataFiles/ibd.intervals.RData")

####### PLOTS ------------------------------------------------------------------
plot.snp.density <- function(snpDensityData, intervals=NULL){

    plotlist <- list()
    chrLens <- trunc(read.table('./dataFiles/chrLengths_mm10.txt')$V2/1e6)

    for (chr in names(snpDensityData)) {
        chrdata <- snpDensityData[[chr]]
        names(chrdata)[1:2] <- c("all", "emp")
        xmax <- chrLens[chr]

        plotlist[[chr]]  <-  ggplot(data=chrdata[[1]],
                                    aes(x=all.Mb, y=Freq, ymin=-100)) +
            geom_bar(stat="identity", fill="#1ab3e3") +     ## medium cyan
            geom_bar(data=chrdata[[2]], aes(x=emp.Mb, y=Freq),
                     stat="identity", fill="#b3e31a") +     ## lime green
            xlab(paste0("Position (Mb) on chromosome ", chr)) +
            ylab("SNP density") +
            scale_x_discrete(limits=(0:xmax), breaks=seq(0, xmax, 10)) +
            theme_bw() +
            theme(axis.title.x = element_text(size=11),
                  axis.title.y = element_text(size=11),
                  axis.text.x = element_text(size=10),
                  axis.text.y = element_text(size=10))

        if (!is.null(intervals)) {
            plotlist[[chr]] <- plotlist[[chr]] +
                annotate("segment", x = intervals[[chr]]$Start,
                         xend = intervals[[chr]]$Stop,
                         y= -40, yend = -40,
                         color="#e31ab3", size=.9, alpha=0.75) ## magenta
        }    # if (!is.NULL(intervals))
    } #  for (chr in seq_along(snpDensityData))
    return(plotlist)
} # END plot.snp.density


### test ----------------------------------------------------------------------
# initialize PDF device
# pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/empSnpDensity.pdf", height=30, width=15, title="Genome-wide SNP Density in the Lg x Sm AIL", colormodel="cmyk")
# # call multiplot and have it plot each chr's data on 1 page
# multiplot(plotlist=plots, cols=1)
# dev.off()

# plots <- plot.snp.density(snpDensityData=snpDensityData, intervals=ibd.intervals)
# pdf(file="./figures/snpDensity_ibdIntervals.pdf", height=30, width=8)
# multiplot(plotlist=plots, cols = 1)
# dev.off()




