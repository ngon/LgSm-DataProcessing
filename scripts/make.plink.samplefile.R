### MAKE SAMPLE FILE FOR PLINK v1.9 ###
### Sample file is in 'oxford format', the format used in IMPUTE2 and related
### tools. See http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html#Sample_File_Format_ for details.

geno.samples<- read.table("./genotyped.samples.txt", as.is=T, header=F)
phenos <- read.table("./dataFiles/phenotypes/allData.txt", sep="\t", as.is=T, header=T)[c(1,43)]
ped <- read.table("./pedigree/pedforQTLRel.txt", sep="\t", as.is=T, header=T)
names(geno.samples) <- "id"
geno.samples <- geno.samples + 0.1

load("./fmiss.RData")
fmiss <- fmiss[c(1,21)]
head(fmiss)


data.tmp <- merge(geno.samples, ail.ped, all.x=TRUE)
data.tmp$sire <- NULL; data.tmp$dam <- NULL
data.tmp$id <- data.tmp$id - 0.1
data.tmp <- merge(data.tmp, phenos, all.x=T)

data <- data.frame(ID_1=data.tmp$id, ID_2=data.tmp$id, missing=fmiss$meanMissing,
                   sex=data.tmp$sex,phenotype=data.tmp$act1.t)
headerline2 <- c(0,0,0,"D","P")
data <- rbind(headerline2, data)
head(data)

write.table(data, file="ailforPlink.sample", sep=" ", row.names=F, col.names=T, quote=F)



