setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles")

# load phenotype data
phenos<- read.table("allData.txt", sep="\t", header=TRUE, as.is=TRUE)

sens1 <- phenos$act4.1 - phenos$act2.1
sens2 <- phenos$act4.2 - phenos$act2.2
sens3 <- phenos$act4.3 - phenos$act2.3
sens4 <- phenos$act4.4 - phenos$act2.4
sens5 <- phenos$act4.5 - phenos$act2.5
sens6 <- phenos$act4.6 - phenos$act2.6

cpp.diff1 <- phenos$cpp8.1 - phenos$cpp1.1
cpp.diff2 <- phenos$cpp8.2 - phenos$cpp1.2
cpp.diff3 <- phenos$cpp8.3 - phenos$cpp1.3
cpp.diff4 <- phenos$cpp8.4 - phenos$cpp1.4
cpp.diff5 <- phenos$cpp8.5 - phenos$cpp1.5
cpp.diff6 <- phenos$cpp8.6 - phenos$cpp1.6

allData <- cbind(phenos, sens1, sens2, sens3, sens4, sens5, sens6,
                 cpp.diff1, cpp.diff2, cpp.diff3, cpp.diff4, cpp.diff5,
                 cpp.diff6)

# write to allData.txt
write.table(allData, "allData.txt", sep="\t", row.names=F, col.names = T)

# write to phenotypes in /dataFiles
phenotypes<- read.table("phenotypes.txt", sep="\t", header=T)
phenotypes<- cbind(phenotypes, sens1, sens2, sens3, sens4, sens5, sens6,cpp.diff1, cpp.diff2, cpp.diff3, cpp.diff4, cpp.diff5, cpp.diff6)
write.table(phenotypes, "phenotypes.txt", sep="\t", row.names=F, col.names=T)

# write to phenotypes and pheno.names in /LgSmDataProcessing
write.table(names(phenotypes), "../pheno.names.txt", sep="\t", row.names=F)
write.table(phenotypes, "../phenotypes.txt", sep="\t", row.names=F, col.names=T)
write.table(phenotypes, "../pheno.noHeader.txt", sep="\t", row.names=F, col.names=F)