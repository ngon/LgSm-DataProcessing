# load positions of empirical snps (a list called empRows)
setwd("/group/palmer-lab/AIL/GBS/genoSummaries")
load("/group/palmer-lab/AIL/GBS/genoSummaries/empRows.Rdata")
source("/group/palmer-lab/AIL/LgSm-DataProcessing/cookieTime/multiplot.R")
library("ggplot2")

### get.empirical.data ----------------------------------------------------
get.empirical.data <- function(file, e=emp.rows, x=maf, y=het.snp){

    print("Getting data on empirical SNPs...")

    chrname             <- substr(basename(file), 1, nchar(basename(file))-16)
    empR                <- e[[chrname]]
    emp.maf             <- data.frame(table(sort(x[seq_along(x) %in% empR])))
    names(emp.maf)[1]   <- "emp.maf"
    emp.het             <- data.frame(table(sort(y[seq_along(y) %in% empR])))
    names(emp.het)[1]   <- "emp.het.snp"

    list(emp.maf=emp.maf, emp.het=emp.het)
}

### maf and het --------------------------------------------------------------
maf.and.het <- function(file, upper=1.2, lower=0.8, emp.rows=NULL) {

    print("Getting genotype data...")
            geno        <- read.table(file, header=F, as.is=T)[-c(1:3)]

    print("Calculating minor allele frequencies...")
            freq        <- cbind( (rowMeans(geno, na.rm = T)/2),
                            1-(rowMeans(geno, na.rm=T)/2)  )
            maf         <- round(apply(freq, 1, min), digits=1)

    print("Calculating average heterozygosity...")
            hetsites    <- geno <= upper & geno >= lower
            het.mouse   <- apply(hetsites, 2, mean)
            het.snp     <- round(apply(hetsites, 1, mean), digits=1)

    if(!is.null(emp.rows)){
        empStats          <- get.empirical.data(file=file, e=emp.rows, x=maf, y=het.snp)

        maf               <- data.frame(table(sort(maf)))
        names(maf)[1]     <- "maf"
        het.snp           <- data.frame(table(sort(het.snp)))
        names(het.snp)[1] <- "het.snp"
        result  <- list(maf, het.snp, het.mouse, empStats["emp.maf"], empStats["emp.het"])

    } else {
        maf               <- data.frame(table(sort(maf)))
        names(maf)[1]     <- "maf"
        het.snp           <- data.frame(table(sort(het.snp)))
        names(het.snp)[1] <- "het.snp"
        result<-  list(maf, het.snp, het.mouse)
    }
    return(result)

}


## get a list of genotype file names
chromosomes <- paste0("chr", 1:19)
filenames <- list()
for (i in chromosomes){
    filenames[i] <- paste0("../dosage/", i, ".filtered.dosage")
}

## get data
snpInfo <- lapply(filenames, maf.and.het, emp.rows=empRows)

## plot data
mafplots <- list()
hetplots <- list()
mouseplots <- list()

data <- snpInfo

for(chr in chromosomes){
    print("Plotting MAF for each chromosome....")
    all <- data.frame(data[[chr]][1])
    names(all)[1:2] <-c("maf", "Freq")
    emp <- data.frame(data[[chr]][4])
    names(emp)[1:2] <-c("maf", "Freq")

    mafplots[[chr]] <- ggplot(data=all, aes(x=maf, y=Freq/sum(Freq))) +
        geom_bar(stat="identity", fill="goldenrod1", color="black") +
        geom_bar(data=emp, aes(x=maf, y=Freq/sum(Freq)),
                 stat="identity", fill="orangered", alpha=0.5, color="black") +
        xlab(paste0("Minor allele frequency on ", chr)) +
        ylab("Proportion of SNPs") +
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=11),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))

    print("Plotting SNP heterozygosity for each chromosome....")
    all <- data.frame(data[[chr]][2])
    names(all)[1:2] <-c("het", "Freq")
    emp <- data.frame(data[[chr]][5])
    names(emp)[1:2] <-c("het", "Freq")

    hetplots[[chr]] <- ggplot(data=all, aes(x=het, y=Freq/sum(Freq))) +
        geom_bar(stat="identity", fill="goldenrod1", color="black") +
        geom_bar(data=emp, aes(x=het, y=Freq/sum(Freq)),
                 stat="identity", fill="orangered", alpha=0.5, color="black") +
        xlab(paste0("Average heterozygosity on ", chr)) +
        ylab("Proportion of SNPs") +
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=11),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))

    print("Plotting mouse heterozygosity for each chromosome....")
    all <- as.data.frame(unlist(data[[chr]][3]))
    names(all)[1] <- "het"

    mouseplots[[chr]] <- ggplot(data=all, aes(x=het)) +
        geom_histogram(fill="goldenrod1", color="black", binwidth=0.1) +
        xlab(paste0("Average heterozygosity on ", chr)) +
        ylab("Proportion of mice") +
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=11),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))

}



pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/mafChromosomes.pdf",
    height=12, width=10, colormodel="cmyk")
multiplot(plotlist=mafplots, cols=4)
dev.off()

pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/hetsnpChromosomes.pdf",
    height=12, width=10, colormodel="cmyk")
multiplot(plotlist=hetplots, cols=4)
dev.off()

pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/hetmouseChromosomes.pdf",
    height=10, width=8, colormodel="cmyk")
multiplot(plotlist=mouseplots, cols=5)
dev.off()
