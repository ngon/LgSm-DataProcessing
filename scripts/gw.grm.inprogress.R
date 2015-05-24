setwd("/group/palmer-lab/AIL/GBS/genoSummaries")
#load("/group/palmer-lab/AIL/GBS/genoSummaries/empRows.Rdata")


gw.grm <- function(file) {
    geno <- as.matrix(read.table(file, header=F)[-c(1:3)])
    meanDosage <- rowMeans(geno, na.rm = T)
    centerGeno <- (scale(t(geno), center=meanDosage, scale=FALSE))
    GRM <- centerGeno%*%(t(centerGeno)) / nrow(centerGeno)
    pdf(file="grm.pdf", width = 7, height=7)
    image(GRM, main="AIL genetic relationship matrix")
}


### plot mafs
chromosomes <- paste0("chr", 1:19)
filenames <- list()
for (i in chromosomes){
    filenames[i] <- paste0("../dosage/", i, ".filtered.dosage")
}





# compute LD
LD <- (t(centerGeno)%*%centerGeno/ncol(centerGeno))**2

# get snp names
snps <- read.table("chr19.txt")[2]
snps <- as.list(substring(snps$V2, first=11))
rownames(LD)[1:500] <- snps
colnames(LD)[1:500] <- snps
