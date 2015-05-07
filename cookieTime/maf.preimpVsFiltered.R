# load positions of empirical snps (a list called empRows)
setwd("/group/palmer-lab/AIL/GBS/genoSummaries")
#source("/group/palmer-lab/AIL/LgSm-DataProcessing/cookieTime/multiplot.R")
#library("ggplot2")

### STEP ONE ----------------------------------------------------------------
# Q: How many SNPs from preimpute.geno files made it into filtered dosage files?
# A: I used filtered.dosage.emp files generated with the script make.emp.files
# to compare the name & number of empirical SNPs in filtered.dosage files to the
# SNPs in the preimpute.geno.

# get preimpute file names and extract SNP names from each file
preimpFiles <- list.files(path="../preimpute/", pattern="*.preimpute.geno")
preimpFiles <- preimpFiles[-20]
get.positions <- function(file) {read.table(file, as.is=T, header=F)[2]}
preSnps <- lapply(file.path("../preimpute/", preimpFiles), get.positions)
save(preSnps, file="./preSnps.Rdata")
preSnps <- data.frame(unlist(preSnps, use.names=FALSE)) # 316,013 SNPs
names(preSnps)[1] <- "snps"
preSnps$snps <- as.character(preSNPs$snps)

# get SNP names from filtered.dosage files (saved as txt file in /genoSummaries)
filteredSnps <- read.table("./snpNames/empirical.snp.list.txt", header=T,
                           sep="\t", as.is=T)
filteredSnps$snps <- as.character(filteredSnps$snps) # 316,013 SNPs

#preImp and filteredSnps are of identical length.
preInFiltered <- which(preimputedList$snps %in% filteredSnps$snps)


### FUNCTION: convert.to.counts ------------------------------------------------
# convert.to.dosage takes the preimpute.geno files, which contain geno probs,
# and convert them to dosage. For each mouse there are 3 cols: probability
# of 0, 1, or 2 reference alleles.
convert.to.dosage <- function (singleMarkerGPs){
    temp = matrix(singleMarkerGPs, nc=3, byrow=T)
    return (as.vector(temp[,3]*2 + temp[,2]))
}


### FUNCTION: maf --------------------------------------------------------------
# maf returns a list for each chromosome: [1] is maf for empirical SNPs in
# filtered.dosage files, [2] is histogram data of maf for plotting. [3-4] are
# the corresponding results for preimpute.geno SNPs.

maf <- function(file) {

    chrname <- substr(basename(file), 1, nchar(basename(file))-20)
    preimp <- paste0("../preimpute/ail.", chrname, ".preimpute.geno")

    print("Getting genotype data...")
    geno    <- read.table(file, header=F, as.is=T)[-c(1:4)]

    genoProbs <- read.table(preimp, header=F, as.is=T)[-c(1:5)]
    pre.geno <- apply(genoProbs, 1, convert.to.dosage)
    pre.geno <- t(pre.geno)

    print("Calculating minor allele frequencies for filtered emp SNPs...")
    freq        <- cbind( (rowMeans(geno, na.rm = T)/2),
                          1-(rowMeans(geno, na.rm=T)/2)  )
    maf         <- round(apply(freq, 1, min), digits=2)
    maf.hist    <- hist(maf, breaks=seq(0,0.5, 0.05), plot=F,
                        include.lowest=T)

    print("Calculating minor allele frequencies for emp SNPs pre-imputation...")
    pre.freq        <- cbind( (rowMeans(pre.geno, na.rm = T)/2),
                              1-(rowMeans(pre.geno, na.rm=T)/2)  )
    pre.maf         <- round(apply(pre.freq, 1, min), digits=2)
    pre.maf.hist    <- hist(pre.maf, breaks=seq(0,0.5, 0.05), plot=F,
                            include.lowest=T)

    result<-  list(maf, maf.hist, pre.maf, pre.maf.hist)
    return(result)

}


### get data-----------------------------------------------------------------
chromosomes <- paste0("chr", 1:19)
empfiles <- list()
for (i in chromosomes){
    empfiles[i] <- paste0("../dosage/", i, ".filtered.dosage.emp")
}

maf.data <- lapply(empfiles, maf)
for(chr in names(maf.data)){
    names(maf.data[[chr]]) <- c("filtered", "filt.hist", "preimp", "pre.hist")
    }
# in most cases the preImp data has a lower #SNPs than the filtered data. I'm
# guessing this is because missing values were removed while calculating MAF,
# and missing SNPs were imputed in filtered.dosage files.




for(chr in chromosomes){
    fileName <- paste0("/group/palmer-lab/AIL/LgSm-DataProcessing/figures/maf.filt",
                       chr, ".pdf")

    pdf(file=fileName,height=4, width=4)

        plot(maf.data[[chr]]$filt.hist, main=paste0("MAF for filtered SNPs on ", chr),
           xlab="MAF", ylab="Number of SNPs", cex.axis=0.9)
    #plot(maf.data[[chr]]$pre.hist, main=paste0("MAF for preimpute SNPs on ", chr),
     #    xlab="MAF", ylab="Number of SNPs", cex.axis=0.9)

    dev.off()
}




plot.maf <- function(all, emp) {
    plot(x=all, xlab=paste0("Minor allele frequency on ", chr),
         ylab="Number of SNPs", main=NULL, cex.lab=0.9, cex.axis=0.8,
         col=rgb(255,193,37, maxColorValue = 255))

    plot(x=emp, col=rgb(255,69,0,128, maxColorValue = 255), add=T)

    legend(x=0, y=max(emp$counts), legend=c("Filtered", "Preimputed"), bty="n",
           cex=0.8, fill=c(rgb(255,193,37, maxColorValue=255),
                           rgb(255,69,0,150, maxColorValue=255)))
    box()

}


# plot maf.snp
pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/mafHistsTest.pdf",
    height=16, width=14.5, colormodel="cmyk")
par(mfrow=c(5,4))
for(chr in chromosomes){
    filt <- maf.data[[chr]]$filt.hist
    pre <- maf.data[[chr]]$pre.hist
    plot.maf(all=filt, emp=pre)
}
dev.off()


















# plot correl between preimp and filtered SNPs
pdf(file="./mafDot_preimpVsFilt.pdf",
    height=12, width=14, colormodel="cmyk")
par(mfrow=c(5,4))

for(chr in chromosomes){
    plot(x=data.frame(filtered[chr]), y=data.frame(preimpute[chr]), type="p", pch=".",
         main="MAF in filtered vs preimputed SNPs",
         ylab="MAF before imputation", xlab="MAF for filtered SNPs")
    }
dev.off()

correlation <- list()
for(chr in chromosomes){
    correlation[chr] <- cor(x=data.frame(filtered[chr]), y=data.frame(preimpute[chr]))
}

