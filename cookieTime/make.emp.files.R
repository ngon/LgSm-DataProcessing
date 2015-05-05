## May 1 2015 ##

## Summary ##
## Instead of running GEMMA on all 4e10+6 SNPs, which can take a long time, I'm
## mapping QTL using empirical SNPs only. I already have a list of the empirical
## SNPs within /GBS/dosage files. Here I'm creating new files for empirical SNPs
## to feed to GEMMA. I want to know the locations and distribution of emp SNPs
## along each chromosome and the read depths associated with them. I will code
## them as LG or SM alleles so we can visualize putatitive recombination points
## and get an idea of how long our haplotypes are. This may help me to identify
## misgenotyped SNPs.

## 1. Make emp-only dosage and snpinfo files for GEMMA ##
setwd("/group/palmer-lab/AIL/GBS/dosage")
load("/group/palmer-lab/AIL/GBS/genoSummaries/empRows.Rdata")

## Function that gets empirical SNPs from complete data set and makes new files.
make.emp.files <- function(dosageFile, empRows){

    geno        <- read.table(dosageFile, header=FALSE, as.is=TRUE)
    chrname     <- substr(dosageFile, 1, nchar(dosageFile)-16)
    emp         <- empRows[[chrname]]
    emp.snps    <- geno[row.names(geno) %in% emp == TRUE,]
    write.table(emp.snps, file=paste0(dosageFile, ".emp"), quote=FALSE, sep="\t",
                col.names=FALSE, row.names=FALSE)

    infoFile    <- paste0(substr(dosageFile, 1, nchar(dosageFile)-6), "snpinfo")
    info        <- read.table(infoFile, header=FALSE, as.is=TRUE)
    emp.info    <- info[row.names(info) %in% emp == TRUE,]
    write.table(emp.info, file=paste0(infoFile, ".emp"), quote=FALSE, sep="\t",
                col.names=FALSE, row.names=FALSE)
}

## Execute the function make.emp.files.
chromosomes <- paste0("chr", 1:19)
filenames <- list()
for (i in chromosomes){ filenames[i]   <- paste0(i, ".filtered.dosage") }
lapply(filenames, make.emp.files, empRows=empRows)









