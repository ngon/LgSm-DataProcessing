
setwd("/home/ngonzales/")

install.packages("QTLRel")
library("QTLRel")

mouselist <- read.table("/group/palmer-lab/AIL/LgSm-DataProcessing/pedigree/mouselist.qtlRel.txt", sep="\t", header=T)

ped <- read.table("/group/palmer-lab/AIL/LgSm-DataProcessing/pedigree/testped.txt", sep="\t", header=T)


ail.idcf <- cic(ped=ped, ids=mouselist$id, df=5, ask=F, verbose=T)

save(list=c("ail.idcf"), file="ail.idcf.RData")
