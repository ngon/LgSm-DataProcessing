### MAKE SAMPLE FILE FOR PLINK v1.9 ###
### Sample file is in 'oxford format', the format used in IMPUTE2 and related
### tools. See http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html#Sample_File_Format_ for details.


# LOAD DATA
setwd("/group/palmer-lab/AIL/LgSm-DataProcessing")
phnames <- read.table(file="./phenotypes.txt",sep="\t", header=T, as.is=T)
pheno<- read.table(file="phenos.allgeno.txt",sep="\t",header=F, as.is=T)[1:30]
names(pheno)[1] <- "id"

covars <- read.table("./covariates.orm.txt", sep="\t", header=T, as.is=T)[1:30]
covars <- data.frame(covars$id, covars$gen, covars$sex)
phnames <- phnames[phnames$id %in% covars$covars.id,]

# i'm using only 1 phenotype column for now. i randomly chose cpp1.1.
# ail mice F50-56 will have values in this col; other mice will have NA.
ails <- cbind(covars, phnames$cpp1.1)
names(ails) <- c("id", "gen", "sex", "cpp1.1")
ails$sex<- as.factor(ails$sex)
levels(ails$sex) <- c(1,2)

## sex and gen data for other AIL mice
othermice<- read.table("C:/Users/Administrator/Desktop/ail39-43sexinfo.txt",
                       sep="\t", header=T, as.is=T, strip.white=T, nrow=840)[1:4]
names(othermice)[4] <- "cpp1.1"
othermice$sex <- as.factor(othermice$sex)
levels(othermice$sex) <- c(2,1) # females=2, males=1

# load genotyped samples
geno.samples <- read.table("./dataFiles/genotyped.samples.txt", sep="\t",
                           as.is=T, header=F)

sam1 <- othermice[othermice$id %in% geno.samples$id,] # 715 samples
sam2 <- ails[ails$id %in% geno.samples$id,] # 1069 individuals

all.not <- geno.samples[!geno.samples$id %in% alldata$id,]
gen<-c(42,42,42,42, 43, 50,50,50,50,50,50,50,50,50,50,50,50,50,51,51,51,51,51,NA,51, 51, 53,53,52,52, 52,52,53,53,53,53,53,54,54,54)
sex<- c(1,1,1,1, 2, 1, 2,2,2,1,1,1,1,1,1,2,2,1,2,2,1,1,2, NA, 2,1, 1,2, 2,1,2,2,2,2,1,1,2,1,1,2)
missingRows <- cbind(as.integer(all.not[1:40]), as.integer(gen), as.integer(sex))


49290**
[41] "54367" "54368" "54371" "54372" # GBS flowcell 27 lib 79, 76, 77, 79
"54386" # flowcell 28 lib 84 - gen 56
"57801" # flowcell 25 lib 65 - allegedly generation 54



test <- merge(geno.samples, alldata, all.x=T)





save(ails, othermice, sampleFile, file="plinkTemp.RData")
sampleFile <- data.frame(ID_1=alldata[1], ID_2=alldata[1], missing=rep(NA, times=1830),
                         sex=alldata[42], gen=alldata[31], cpp1.1=alldata$V22)
names(sampleFile)[1:6] <- c("ID_1", "ID_2", "missing", "sex", "gen", "cpp1.1")
sampleFile[1,] <- c(0, 0, 0, "D", "D","P"), sampleFile)

# save file
write.table(sampleFile, file="plink.ailSample", sep=" ", col.names=T, row.names=F,
            quote=F)



