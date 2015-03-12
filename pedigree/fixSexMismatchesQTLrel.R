setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/pedigree")

library("QTLRel")

# read in original phenotype file to get the correct sex information from test mice (e.g. the mice in 'mouse list')
pheno <- read.table("../dataFiles//allData.txt", sep="\t", header=T)
pheno <- cbind(pheno[1:2],pheno[5:6], pheno[13])
levels(pheno$sex)[levels(pheno$sex)=="M"] <- 1
levels(pheno$sex)[levels(pheno$sex)=="F"] <- 2

# read in list of test mice and reverse their sex information to match coding in pheno (M=1, F=2)
mouselist <- read.table("qtlReltestmiceGen.csv", sep=",", header=T)
mouselist$sex <- as.factor(mouselist$sex)
levels(mouselist$sex)[levels(mouselist$sex) == "1"] <- "F"
levels(mouselist$sex)[levels(mouselist$sex) == "2"] <- "1"
levels(mouselist$sex)[levels(mouselist$sex) == "F"] <- "2"

# read in pedigree
origped <- read.table("qtlRelpedGen.csv",sep=",", header=T )
head(origped)



# try to run qtlrel's cic command - it will fail but it will print out ids with sex mismatches. 
idcfs <- cic(ped=origped, ids=mouselist$id, df=5, ask=T, verbose=T)

# generation 43 and 48 has the sires and dams mixed up in its AILineInfo file; fix this by switching the columns for this generation and returning it to the pedigree data frame. 

gen43 <- origped[origped$generation == 43,]
gen48 <- origped[origped$generation == 48,]

gen43 <- cbind(gen43[1], gen43[3], gen43[2], gen43[4:5])
gen48 <- cbind(gen48[1], gen48[3], gen48[2], gen48[4:5])
names(gen43)[1:5]<-names(origped[1:5])
names(gen48)[1:5]<-names(origped[1:5])

head(gen43); tail(gen43)
head(gen48); tail(gen48)

ped <- rbind(origped[1:6389,], gen43, origped[6574:7173,], gen48, origped[7439:10177,])



# test cic command again 
idcfs <- cic(ped=ped, ids=mouselist$id, df=5, ask=T, verbose=T)

# 46950 is a dam but is listed as a sire in the origped/ped for individual 48705 [3676,]. The sire is 46359. 
ped[ped$sire == 46950,]
ped[3676,2]<-46359
ped[3676,3]<-46950

ped.save <- ped
# generations 46 and 51 have mix-ups as well - but not in all cases. therefore i can't just switch the dam and sire columns in generations 47 and 52. instead i'm going to read in a list of the mismatched animals and 

# get list of female mice from gen 46 and 51 that are incorrectly listed as sires in generations 47 and 52.
mism <- read.table("sex.mismatches.qtlrel.txt", header=F)
mism <- mism[5:7]
names(mism)<-c("sex", "generation", "id")

# match mism$id to sire$ids in ped and switch dam/sire cols
fixthese <- ped[ped$sire %in% mism$id,]
fixthese <- cbind(fixthese[1], fixthese[3], fixthese[2], fixthese[4:5])
names(fixthese)[2:3]<-c("true.sire", "true.dam")
head(fixthese)

# match and replace rows in ped with corresponding rows in fixthese
library("plyr")

#indexes <- match(fixthese$id, ped$id)
ped$id <- as.character(ped$id)
fixthese$id <- as.character(fixthese$id)

# i merged ped and fixthese so that the true.dam and true.sire were listed in extra columns in testped. all non-matches had an NA in those columns. I then opened the file in excel and corrected the dam/sire info by replacement. 
testped <- merge(ped, fixthese, by="id", all.y=TRUE, all.x=TRUE)
testped <- testped[-c(8:9)]
write.table(testped, "testped.txt", sep="\t", col.names=T, row.names=F)
write.table(mouselist, "mouselist.qtlRel.txt", sep="\t", col.names=T, row.names=F)

# try running cic command on fixedped
testped <- read.table("testped.txt", sep="\t", header=T)


# test cic command again 
idcfs <- cic(ped=testped, ids=mice$id, df=5, ask=T, verbose=T)

# success! disk space needed = 83.2 Gb
# however there are still errors... program says certain IDs are out of bounds, but they shouldn't be.

# import inbred ped fragment from peter's lgsmfear githib repo
# except wtf peter... the sire and dams are identical. this is an error. 
inbred.ped <- read.table("pedigreeInbred.txt", sep="\t", header=T)
inbred.ped <- cbind(inbred.ped[1], inbred.ped[4:5], inbred.ped[3], inbred.ped[2])
levels(inbred.ped$sex)[levels(inbred.ped$sex)=="M"] <- "1"
levels(inbred.ped$sex)[levels(inbred.ped$sex)=="F"] <- "2"
write.table(inbred.ped, "inbredpedNMG.txt", sep="\t", row.names=F)

testped <- read.table("testped.txt", sep="\t", header=T)
testped <- testped[1:5]
testped <- testped[3:10177,]
testped$id <- as.factor(testped$id)
testped$sire <- as.factor(testped$sire)
testped$dam <- as.factor(testped$dam)
testped$sex <- as.factor(testped$sex)
testped$generation <- as.factor(testped$generation)
write.table(testped, "pedforQtlRel.txt", sep="\t", row.names=F)

ped <- read.table("pedforQtlRel.txt", sep="\t", header=T)
mice$id <- as.factor(mice$id)
write.table(mice, "mice.txt", sep="\t", row.names=F)

ail.kinship <- kinship(ped, mice$id)

# test cic command again 
idcfs <- cic(ped=ped, ids=mice$id, df=5, ask=T, verbose=T)
