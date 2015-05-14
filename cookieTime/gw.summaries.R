# this code takes approximately 9 hours to run on CRI when 150 GB are allocated
#setwd("/group/palmer-lab/AIL/GBS/genoSummaries")
#load("/group/palmer-lab/AIL/GBS/genoSummaries/empRows.Rdata")
library("ggplot2")

### genome wide maf --------------------------------------------------------------
### all snps
gw.maf <- function(file) {
    geno <- read.table(file, header=F, as.is=T)[-c(1:3)]
    freq <- cbind((rowMeans(geno, na.rm = T)/2),
                   1-(rowMeans(geno, na.rm=T)/2))
    maf <- round(apply(freq, 1, min), digits=2)
    return(maf)
}

### empirical snps
# gw.emaf <- function(file, e=emp.rows){
#
#     geno <- read.table(file, header=F, as.is=T)[-c(1:3)]
#     freq <- cbind((rowMeans(geno, na.rm = T)/2),
#                   1-(rowMeans(geno, na.rm=T)/2))
#     maf <- round(apply(freq, 1, min), digits=2)
#     chrname <- substr(basename(file), 1, nchar(basename(file))-16)
#     #chrname <- substr(basename(file), 1, nchar(basename(file))-4)
#     empR <- e[[chrname]]
#     emp.maf <- sort(maf[seq_along(maf) %in% empR])
#     return(emp.maf)
# }

### plot mafs
chromosomes <- paste0("chr", 1:19)
filenames <- list()
for (i in chromosomes){
    filenames[i] <- paste0("../dosage/onlyEmpirical", i, ".filtered.dosage")
}

all_maf <- data.frame(table(sort(unlist(lapply(filenames, gw.maf)))))
all_maf$Var1 <- as.numeric(as.character(all_maf$Var1))
# emp_maf <- data.frame(table(sort(unlist(lapply(filenames, gw.emaf, e=empRows)))))
# emp_maf$Var1 <- as.numeric(as.character(emp_maf$Var1))

pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/gwMAF.onlyEmp.pdf",
    height=4, width=4, colormodel="cmyk")

ggplot(data=all_maf, aes(x=Var1, y=Freq/sum(Freq))) +
    geom_bar(stat="identity", fill="goldenrod1", color="black") +
    #geom_bar(data=emp_maf, aes(x=Var1, y=Freq/sum(Freq)), stat="identity",
     #                          fill="grey25", alpha=0.5, color="black") +
    scale_x_continuous(breaks=seq(0,0.5, 0.05))+
    xlab("Minor allele frequency") +
    ylab("Proportion of SNPs") +
    theme_bw() +
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10))

dev.off()




### genome wide het --------------------------------------------------------------
gw.het <- function(file, upper=1.2, lower=0.8) {
    geno <- read.table(file, header=F, as.is=T)[-c(1:3)]
    hetsites <- geno <= upper & geno >= lower
    het.snp <- round(apply(hetsites, 1, mean), digits=2)
    return(het.snp)
    }

# gw.ehet <- function(file, upper=1.2, lower=0.8, e=emp.rows){
#     geno <- read.table(file, header=F, as.is=T)[-c(1:3)]
#     hetsites <- geno <= upper & geno >= lower
#     het.snp <- round(apply(hetsites, 1, mean), digits=2)
#     chrname <- substr(basename(file), 1, nchar(basename(file))-16)
#     #chrname <- substr(basename(file), 1, nchar(basename(file))-4)
#     empR <- e[[chrname]]
#     emp.het <- sort(het.snp[seq_along(het.snp) %in% empR])
#     return(emp.het)
# }



all_het <- data.frame(table(sort(unlist(lapply(filenames, gw.het)))))
all_het$Var1 <- as.numeric(as.character(all_het$Var1))
# emp_het <- data.frame(table(sort(unlist(lapply(filenames, gw.ehet, e=empRows)))))
# emp_het$Var1 <- as.numeric(as.character(emp_het$Var1))

pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/gwHET.onlyEmp.pdf",
    height=4, width=4, colormodel="cmyk")

ggplot(data=all_het, aes(x=Var1, y=Freq/sum(Freq))) +
    geom_bar(stat="identity", fill="darkolivegreen3", color="black") +
#     geom_bar(data=emp_het, aes(x=Var1, y=Freq/sum(Freq)), stat="identity",
#              fill="grey25", alpha=0.5, color="black") +
    scale_x_continuous(breaks=seq(0,0.5, 0.05))+
    xlab("Average SNP heterozygosity") +
    ylab("Proportion of SNPs") +
    theme_bw() +
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10))

dev.off()







