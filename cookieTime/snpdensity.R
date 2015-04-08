

###################################################
setwd("/group/palmer-lab/AIL/GBS/genoSummaries/")
library("ggplot2")
source("/group/palmer-lab/AIL/LgSm-DataProcessing/cookieTime/multiplot.R")

# load list of empirical SNP positions (emp)
# and chromosome lengths (chrlens)

load("./empSnpPositions.RData")

snp.density <- function(filename, empRows=FALSE) {
    chrPos <- trunc((read.table(filename, sep="\t", header=F)[2])/1e6)
    allDensity <- data.frame(table(sort(chrPos$V2)))
    names(allDensity)[1] <- "all.Mb"

    if(empRows==TRUE){
        is.emp <- which((read.table(filename, sep="\t", header=F, as.is=T)[2])$V2 %in% emp$ps)
        empDensity <- data.frame(table(sort(chrPos[row.names(chrPos) %in% is.emp,])))
        names(empDensity)[1] <- "emp.Mb"
        list(allDensity, empDensity)
    } else {
        list(allDensity)
    }

}

####### COMMANDS  --------------------------------------------------------------
chromosomes <- paste0("chr", 1:19)

filenames <- list()
for (i in chromosomes){
    filenames[i] <- paste0("../dosage/", i, ".filtered.snpinfo")
}

snpsPerMb <- list()
for(file in filenames) {
    snpsPerMb[[file]] <- snp.density(file, empRows=TRUE)
}
names(snpsPerMb)[1:19] <- paste0("chr", 1:19)
save(snpsPerMb, file="./snpsPerMb.RData")

####### PLOTS ------------------------------------------------------------------
plotlist <- list()

for (chr in seq_along((snpsPerMb))) {
    chrdata <- snpsPerMb[[chr]]
    names(chrdata)[1:2] <- c("all", "emp")
    xmax <- chrlens[chr,]

    plotlist[[chr]] <- ggplot(data=chrdata[[1]], aes(x=all.Mb, y=Freq)) +
        geom_bar(stat="identity", fill="goldenrod1") +
        geom_bar(data=chrdata[[2]], aes(x=emp.Mb, y=Freq),
                 stat="identity", fill="steelblue3") +
        xlab(paste0("Position (Mb) on chromosome ", chr)) +
        ylab("SNP density") +
        scale_x_discrete(limits=(0:xmax), breaks=seq(0, xmax, 10)) +
        #scale_y_discrete(limits=(0:400), breaks=seq(0,400, 50))+
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=11),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))
}

pdf(file="./snpDensityPlots.pdf", height=14, width=10, title="Genome-wide SNP Density in the Lg x Sm AIL", colormodel="cmyk")
multiplot(plotlist=plotlist, cols=2)
dev.off()






