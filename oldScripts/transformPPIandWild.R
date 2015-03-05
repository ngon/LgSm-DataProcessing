setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles")

# load phenotype data
data<- read.table("allData.txt", sep="\t", header=TRUE, as.is=TRUE)



### Transform ppi phenotypes ###
################################

# ppi phenos are proportions. transform to the log-odds scale using the following functions:

logit10 <- function(x)
        log10((x + .Machine$double.eps)/(1 - (x + .Machine$double.eps)))
# .Machine$double.eps = smallest positive float x such that 1 + x != 1. Here it is 2.220446e-16. 

project.onto.interval <- function(x, a, b)
        pmin(b, pmax(a,x))

transform.proportion <- function(x) 
        logit10((project.onto.interval(x, 0.01, 0.99))

# transform the phenotypes in allData file 
# values are expresed as percentages; change to proportions first.
ppi3.logit <- transform.proportion(data$ppi3/100)
ppi6.logit <- transform.proportion(data$ppi6/100)
ppi12.logit <- transform.proportion(data$ppi12/100)


### Transform wild phenotype ###
################################

wild.binary <- data$wild
wild.binary <- replace(wild.binary, wild.binary > 0, 1)


### Save all files           ###
################################ 
# write to allData file in /dataFiles
allData <- cbind(data, ppi3.logit, ppi6.logit, ppi12.logit, wild.binary)
write.table(allData, "allData.txt", sep="\t", row.names=F, col.names = T)

# write to phenotypes in /dataFiles
phenotypes<- read.table("phenotypes.txt", sep="\t", header=T)
phenotypes<- cbind(phenotypes, ppi3.logit, ppi6.logit, ppi12.logit, wild.binary)
write.table(phenotypes, "phenotypes.txt", sep="\t", row.names=F, col.names=T)

# write to phenotypes and pheno.names in /LgSmDataProcessing
write.table(names(phenotypes), "../pheno.names.txt", sep="\t", row.names=F)
write.table(phenotypes, "../phenotypes.txt", sep="\t", row.names=F, col.names=T)
write.table(phenotypes, "../pheno.noHeader.txt", sep="\t", row.names=F, col.names=F)


