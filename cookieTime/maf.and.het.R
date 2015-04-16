# load positions of empirical snps (a list called empRows)
setwd("/group/palmer-lab/AIL/GBS/genoSummaries")
load("/group/palmer-lab/AIL/GBS/genoSummaries/empRows.Rdata")
#source("/group/palmer-lab/AIL/LgSm-DataProcessing/cookieTime/multiplot.R")
#library("ggplot2")

### get.empirical.data ----------------------------------------------------
get.empirical.data <- function(file, e=emp.rows, x=maf, y=het.snp){

    print("Getting data on empirical SNPs...")

    chrname             <- substr(basename(file), 1, nchar(basename(file))-16)
    empR                <- e[[chrname]]
    emp.maf             <- sort(x[seq_along(x) %in% empR])
    empmaf.hist         <- hist(emp.maf[1], breaks=seq(0,0.5, 0.05), plot=F,
                        include.lowest=T)
    emp.het             <- sort(y[seq_along(y) %in% empR])
    emphet.hist         <- hist(emp.het[1], breaks=seq(0,0.5, 0.05), plot=F,
                        include.lowest=T)

    list(empmaf.hist=empmaf,hist, emphet.hist=emphet.hist)
}

### maf and het --------------------------------------------------------------
maf.and.het <- function(file, upper=1.2, lower=0.8, emp.rows=NULL) {

    print("Getting genotype data...")
            geno        <- read.table(file, header=F, as.is=T)[-c(1:3)]

    print("Calculating minor allele frequencies...")
            freq        <- cbind( (rowMeans(geno, na.rm = T)/2),
                            1-(rowMeans(geno, na.rm=T)/2)  )
            maf         <- round(apply(freq, 1, min), digits=2)
            maf.hist    <- hist(maf, breaks=seq(0,0.5, 0.05), plot=F,
                            include.lowest=T)

    print("Calculating average heterozygosity...")
            hetsites    <- geno <= upper & geno >= lower
            het.mouse   <- apply(hetsites, 2, mean)
            het.snp     <- round(apply(hetsites, 1, mean), digits=2)
            hetsnp.hist <- hist(het.snp, breaks=seq(0,0.5, 0.05), plot=F,
                            include.lowest=T)

    if(!is.null(emp.rows)){
        empStats <- get.empirical.data(file=file, e=emp.rows, x=maf, y=het.snp)
        result<- list(maf.hist, hetsnp.hist, het.mouse,
                      empStats["empmaf.hist"], empStats["emphet.hist"])
    } else {
        result<-  list(maf, het.snp, het.mouse)
    }
    return(result)

}

### plotting functions ------------------------------------------------------
plot.maf.snp <- function(all, emp) {
    plot(x=all[[1]], xlab=paste0("Minor allele frequency on ", chr),
         ylab="Number of SNPs", main=NULL, cex.lab=1, cex.axis=.9,
         col=rgb(255,193,37, maxColorValue = 255))

    plot(x=emp[[1]]$empmaf.hist, col=rgb(255,69,0,128, maxColorValue = 255), add=T)

    legend(x=0, y=max(all[[1]]$counts), legend=c("All", "GBS"), bty="n",
           cex=0.8, fill=c(rgb(255,193,37, maxColorValue=255),
                           rgb(255,69,0,150, maxColorValue=255)))
    box()
}


plot.het.snp <- function(all, emp) {
    plot(x=all[[1]], xlab=paste0("Average heterozygosity on ", chr),
         ylab="Number of SNPs", main=NULL, cex.lab=1, cex.axis=.9,
         col=rgb(255,193,37, maxColorValue = 255))

    plot(x=emp[[1]]$emphet.hist, col=rgb(255,69,0,128, maxColorValue = 255), add=T)

    legend(x=0, y=max(all[[1]]$counts), legend=c("All", "GBS"), bty="n",
           cex=0.8, fill=c(rgb(255,193,37, maxColorValue=255),
                           rgb(255,69,0,150, maxColorValue=255)))
    box()
}

plot.het.mouse <- function(mice) {
    hist(x=mice[[1]], xlab=paste0("Average heterozygosity on ", chr),
         ylab="Number of mice", main=NULL, cex.lab=1, cex.axis=.9,
         col=rgb(255,193,37, maxColorValue = 255))
    box()
}


### get data-----------------------------------------------------------------
chromosomes <- paste0("chr", 1:19)
filenames <- list()
for (i in chromosomes){
    filenames[i] <- paste0("../dosage/", i, ".filtered.dosage")
}
data <- lapply(filenames, maf.and.het, emp.rows=empRows)


## plot maf.snp
pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/maf_chr.pdf",
    height=12, width=14, colormodel="cmyk")
par(mfrow=c(5,4))
for(chr in chromosomes){
    all <- data[[chr]][1]
    emp <- data[[chr]][4]
    plot.maf.snp(all, emp)
}
dev.off()

## plot het.snp
pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/snphet_chr.pdf",
    height=12, width=14, colormodel="cmyk")
par(mfrow=c(5,4))
for(chr in chromosomes){
    all <- data[[chr]][2]
    emp <- data[[chr]][5]
    plot.het.snp(all, emp)
}
dev.off()

## plot het mouse
pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/mousehet_chr.pdf",
    height=15, width=12, colormodel="cmyk")
par(mfrow=c(5,4))
for(chr in chromosomes){
    mice <- data[[chr]][3]
    plot.het.mouse(mice)
}
dev.off()








