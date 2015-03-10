
setwd("/home/ngonzales/")

library("QTLRel")

mouselist <- read.table("mouselist.qtlRel.txt", sep="\t", header=T)

pedigree <- read.table("testped.txt", sep="\t", header=T)


ail.idcf <- cic(ped=ped, ids=mouselist$id, df=5, ask=T, verbose=T)

save.image()

