### 15 May 2015

### PURPOSE: Add cM positions to univariate GEMMA results and make manhattan
### plots with cM on the x-axis. Also make zoom plots of top associations for
### each trait.

### FILE PATHS
## 1. chrN.bp.cM.txt (\t) files: chr | bp position | cM position; or RData file
##      /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/centimorgans/
##      /group/palmer-lab/AIL/LgSm-DataProcessing/dataFiles/snpChrPosCm.RData
## 2. trait.chrN.assoc.txt
##      /group/palmer-lab/AIL/qtlmapping/output/
## 3. AIL.traits.RData : list with trait names and descriptions
##      /group/palmer-lab/AIL/LgSm-DataProcessing/dataFiles/AIL.traits.RData

### I. ADD CENTIMORGAN COLUMN TO ASSOC FILES -----------------------------------

### get trait names on CRI or from Natalia's local directory
load('/group/palmer-lab/AIL/LgSm-DataProcessing/dataFiles/AIL.traits.RData')
#load("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles/AIL.traits.RData")

### get position and cM information from CRI or from Natalia's local directory
### this gives you a list 'snpinfo' of data frames w/ positions for each chr
load('/group/palmer-lab/AIL/LgSm-DataProcessing/dataFiles/snpChrPosCm.RData')
#load("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles/snpChrPosCm.RData")


append.cM <- function(trait, snpinfo=snpinfo){

=======
    for(trait in traits){
    chromosomes <- paste0("chr", 1:19)

    for (chr in chromosomes){
    assoc      <- read.table(file=paste0("/group/palmer-lab/AIL/qtlmapping/output/",
                                          trait, ".", chr, ".assoc.txt"),
                              header=TRUE, as.is=TRUE, sep="\t")
    whichChr    <- snpinfo[[chr]]
    whichSnps   <- merge(assoc,whichChr, by.x="ps", by.y="bp", all.x=TRUE)
    whichSnps   <- data.frame(whichSnps[c(2,3,1, 4:9,11)])

    write.table(whichSnps, file=paste0("/group/palmer-lab/AIL/qtlmapping/output/cM.",
                                       trait, ".", chr, ".assoc.txt"),
                col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
    }
    }



