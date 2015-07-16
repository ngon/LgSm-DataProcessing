###################################################
library("ggplot2")
source("/group/palmer-lab/AIL/LgSm-DataProcessing/scripts/multiplot.R")
# load list of empirical SNP positions (emp) and chromosome lengths (
load("./empSnpsPosByChr.RData") # object: empSnpsByChr

}

snp.density <- function(chr, empRows=FALSE, empPositions=NULL) {

    filename <- paste0("/group/palmer-lab/AIL/GBS/dosage/", chr, ".filtered.snpinfo")
    chrPos <- trunc((read.table(filename, sep="\t", header=F)[2])/1e6)
    allDensity <- data.frame(table(sort(chrPos$V2)))
    names(allDensity)[1] <- "all.Mb"

    if(empRows==TRUE){
        is.emp <- which((read.table(filename, sep="\t", header=F)[2])$V2 %in% empPositions[[1]])
        empDensity <- data.frame(table(sort(chrPos[row.names(chrPos) %in% is.emp,])))
        names(empDensity)[1] <- "emp.Mb"
        list(allDensity, empDensity)
    } else {
        list(allDensity)
    }

}

####### COMMANDS  --------------------------------------------------------------
chromosomes <- paste0("chr", 1:19)
# list filenames of SNPs from full data set
# filenames <- list()
# for (i in chromosomes){
#     filenames[i] <- paste0("/group/palmer-lab/AIL/GBS/dosage/", i, ".filtered.snpinfo")
# }
snpDensityData <- list()
for(c in chromosomes) {
    snpDensityData[[c]] <- snp.density(chr=c, empRows=TRUE, empPositions = empSnpsByChr[[c]])
}
names(snpDensityData)[1:19] <- paste0("chr", 1:19)
save(snpDensityData, file="/group/palmer-lab/AIL/qtlmapping/snpDensity.allData.RData")


####### PLOTS ------------------------------------------------------------------
chrLens <- trunc(read.table('/group/palmer-lab/AIL/LgSm-DataProcessing/dataFiles/chrLengths_mm10.txt')$V2/1e6)
plotlist <- list()

for (chr in seq_along(snpDensityData)) {
    chrdata <- snpDensityData[[chr]]
    names(chrdata)[1:2] <- c("all", "emp")
   # names(chrdata)[1] <- "all"
    xmax <- chrLens[chr]

    plotlist[[chr]] <- ggplot(data=chrdata[[1]], aes(x=all.Mb, y=Freq)) +
        geom_bar(stat="identity", fill="goldenrod1") +
        geom_bar(data=chrdata[[2]], aes(x=emp.Mb, y=Freq),
                 stat="identity", fill="steelblue3") +
        xlab(paste0("Position (Mb) on chromosome ", chr)) +
        ylab("SNP density") +
        scale_x_discrete(limits=(0:xmax), breaks=seq(0, xmax, 10)) +
        theme_bw() +
        theme(axis.title.x=element_text(size=11),
              axis.title.y=element_text(size=11),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))
}

pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/empSnpDensity.pdf", height=14, width=12, title="Genome-wide SNP Density in the Lg x Sm AIL", colormodel="cmyk")
multiplot(plotlist=plotlist, cols=2)
dev.off()






