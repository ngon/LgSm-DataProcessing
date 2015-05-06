# load positions of empirical snps (a list called empRows)
setwd("/group/palmer-lab/AIL/GBS/genoSummaries")
#source("/group/palmer-lab/AIL/LgSm-DataProcessing/cookieTime/multiplot.R")
#library("ggplot2")

### convert.to.counts ---------------------------------------------------------
convert.to.dosage <- function (singleMarkerGPs){

    temp = matrix(singleMarkerGPs, nc=3, byrow=T)
    return (as.vector(temp[,3]*2 + temp[,2]))
}


### maf  --------------------------------------------------------------
# maf.preimpVsFiltered <- function(file, preimp) {
#
#     print("Getting genotype data...")
#     geno        <- read.table(file, header=F, as.is=T)[-c(1:4)]
#
#     genoProbs <- read.table(preimp, header=F, as.is=T)[-c(1:5)]
#     pre.geno <- apply(genoProbs, 1, convert.to.counts)
#     pre.geno <- t(pre.geno)
#
#
#     print("Calculating minor allele frequencies for filtered emp SNPs...")
#     freq        <- cbind( (rowMeans(geno, na.rm = T)/2),
#                           1-(rowMeans(geno, na.rm=T)/2)  )
#     maf         <- round(apply(freq, 1, min), digits=2)
#     maf.hist    <- hist(maf, breaks=seq(0,0.5, 0.05), plot=F,
#                         include.lowest=T)
#
#     print("Calculating minor allele frequencies for emp SNPs pre-imputation...")
#     pre.freq        <- cbind( (rowMeans(pre.geno, na.rm = T)/2),
#                           1-(rowMeans(pre.geno, na.rm=T)/2)  )
#     pre.maf         <- round(apply(pre.freq, 1, min), digits=2)
#     pre.maf.hist    <- hist(pre.maf, breaks=seq(0,0.5, 0.05), plot=F,
#                         include.lowest=T)
#
#     result<-  list(maf, maf.hist, pre.maf, pre.maf.hist)
#     return(result)
#
# }

### plotting functions ------------------------------------------------------
plot.maf.snp <- function(all, emp) {
    plot(x=all[[1]], xlab=paste0("Minor allele frequency on ", chr),
         ylab="Number of SNPs", main=NULL, cex.lab=1, cex.axis=.9,
         col=rgb(255,193,37, maxColorValue = 255))

    legend(x=0, y=max(all[[1]]$counts), legend=c("Preimpute", "Filtered"), bty="n",
           cex=0.8, fill=c(rgb(255,193,37, maxColorValue=255),
                           rgb(255,69,0,150, maxColorValue=255)))

    plot(x=emp[[1]]$empmaf.hist, col=rgb(255,69,0,128, maxColorValue = 255), add=TRUE)
    box()

}


# load list of SNP positions
load("/group/palmer-lab/AIL/GBS/genoSummaries/snpNames/empSnpsPosByChr.Rdata")
empSnps <- data.frame(unlist(empSnps, use.names=FALSE))
names(empSnps)[1] <- "snps"

# get preimpute filenames except for chrX (idx 20), which is empty
preimpFiles <- list.files(path="../preimpute/", pattern="*.preimpute.geno")
preimpFiles <- preimpFiles[-20]

# extract SNP positions from each file
get.positions <- function(file) {read.table(file, as.is=T, header=F)[3]}
preSnps <- lapply(file.path("../preimpute/", preimpFiles), get.positions)
preSnps <- data.frame(unlist(preSnps, use.names=FALSE))
#preSnps <- do.call(what=rbind.data.frame, args=preSnps) # 316,013 SNPs
names(preSnps)[1] <- "snps"

# get separate lists of empirical and imputed SNPs
empRowsInPreimpute <- which(empSnps$snps %in% preSnps$snps) # 312516
save(preSnps,empRowsInPreimpute, file="./preImputedSnpInfo.RData")



maf.dotplot <- function(file, preimp, emprows) {

    print("Getting genotype data...")
    geno    <- read.table(file, header=F, as.is=T)[-c(1:4)]
    geno    <- geno[which(row.names(geno) %in% emprows),]

    genoProbs <- read.table(preimp, header=F, as.is=T)[-c(1:5)]
    pre.geno <- apply(genoProbs, 1, convert.to.counts)
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
preimpfiles <- list()

for (i in chromosomes){
    empfiles[i] <- paste0("../dosage/", i, ".filtered.dosage.emp")
    preimpfiles[i] <- paste0("../preimpute/ail.", i, ".preimpute.geno")
}

data<- mapply(maf.dotplot, empfiles, preimpfiles, MoreArgs=list(emprows=empRowsInPreimpute))



## plot maf.snp
pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/mafHist_preimpVsFilt.pdf",
    height=12, width=14, colormodel="cmyk")
par(mfrow=c(5,4))
for(chr in chromosomes){
    all <- data[[chr]][2]
    emp <- data[[chr]][4]
    plot.maf.snp(all, emp)
}
dev.off()

# plot correl between preimp and filtered SNPs
pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/mafDot_preimpVsFilt.pdf",
    height=12, width=14, colormodel="cmyk")
par(mfrow=c(5,4))
for(chr in chromosomes){
    plot(x=data[[chr]][1], y=data[[chr]][3], pch=".",
            main=paste("MAF in filtered vs preimputed SNPs\n", cor(data[[chr]][1],
                                                                   data[[chr]][3]), sep=" "),
         ylab="MAF before imputation", xlab="MAF for filtered SNPs")
}
dev.off()

