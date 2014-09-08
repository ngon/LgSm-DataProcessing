# Format pheno files for meeting with Shyam on 8/15/2014
library(xlsx)
setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing")
list.files()

### Read data

# cpp data
cppData<- read.table("./cppPhenos.txt", sep="\t", header=TRUE)

# physiological phenos
other <- read.table("./otherPhenos.txt", sep="\t", na.strings=NA, header=TRUE)
other <- other[do.call(order, other),]


# read ppi data and info sheet separately and merge files
ppiInfo <- read.table(file="ppi.info.csv", sep=",", header=TRUE, na.strings="NA")
ppiInfo <- ppiInfo[do.call(order, ppiInfo),]

ppiData <- read.table(file="ppiJul2014.csv", sep=",", header=TRUE, na.strings="NA")
ppiData <- ppiData[do.call(order, ppiData),]
ppiData <- ppiData[,c(1:2, 81:87)]
names(ppiData)[1:6] <- c("id", "ppi.box", "ppi3", "ppi6", "ppi12", "avg.ppi")

ppiMerge <- merge(ppiInfo, ppiData, by="id", all=TRUE)
which(duplicated(ppiMerge$id)) # there is one duplicated value
ppiMerge[1002:1005,] # 53420 - remove the second entry 
ppiMerge <-ppiMerge[-c(1005),]

# merge ppiMerged data with cpp and other pheno data

ppiOther <- merge(ppiMerge, other, by="id", all=TRUE)
which(duplicated(ppiOther$id)) # 52000
ppiOther[810:813,] # looks like this ID was mixed up in the glucose file
                   # DOUBLE CHECK THIS
# row 812 is the messed up one. it looks like the glucose levels for 52000
# and 51905 were mixed up in this case. I am deleting row 812. 

allData <- merge(cppData, ppiOther, by="id", all=TRUE)
write.table(allData, "./mergedData.txt", sep="\t", row.names=FALSE, quote=FALSE)

# remove incorrect entry for animal 52000 (row 812)
allData <- allData[-c(812),]
# 51329 is also duplicated. merge values and discard NAs.
gluData51329 <- allData[595, 137:144]
allData[594, 137:144] <- gluData51329
allData <- allData[-c(595),]

# sex, gen, fam/sire/dam, and cage variables x and y do not match in all cases
# this is because in one dfs the variable was missing, and in the other it was not
# looking at sex, it appears that this is because the animals were "pheno errors"
# and will not be used to map QTLs. 

blankValues <- which(is.na(allData$sex.x))

allData[allData$id == "46003",] # remove this mouse (row 25)
allData[allData$id == "46368",] # remove
allData[allData$id == "46445",] # remove
allData[allData$id == "49003",] # remove
allData[allData$id == "49290",] # remove
allData[allData$id == "50112",] # remove
allData[allData$id == "50714",] # remove

allData <- allData[-c(blankValues),] # removing entries above
allData[allData$id == "50716",] # fix (row 497) - i'm doing this in excel

# read in new, corrected table
allData<- read.table("./mergedData.txt", sep="\t", header=T, na.strings="NA")

# correct fam mismatches
allData[allData$fam.y == "BrF52-55", 1:13] 
allData[c(539,542), c(11,12)] <- allData[c(539,542), c(5,6)]

# now all duplicated columns are identical! keep one and rename:
allData <- allData[-c(8:13)]
names(allData)[2:7] <- c("gen", "cage", "fam", "dam", "sire", "sex")

# make separate coat color categories
unique(allData$cc)
allData[allData$cc == "AGR/L",] # 587
allData[587, 8] <- "AGRL"
allData[allData$cc == "ARG",]
allData[535, 8] <- "AGR"
allData[95,1:10]
allData[95,8] <- "BL"
allData[93,1:8]
allData[93,8] <- "AGR"
allData$cc <- droplevels(allData$cc)

colors <- allData$cc
levels(colors)[levels(colors)=="AGL"] <- "A"
levels(colors)[levels(colors)=="AGR"] <- "A"
levels(colors)[levels(colors)=="AGRL"] <- "A"
levels(colors)[levels(colors)=="WL"] <- "W"
levels(colors)[levels(colors)=="WR"] <- "W"
levels(colors)[levels(colors)=="WRL"] <- "W"
levels(colors)[levels(colors)=="BL"] <- "B"
levels(colors)[levels(colors)=="BR"] <- "B"
levels(colors)[levels(colors)=="BRL"] <- "B"

allData$agouti <- colors

write.table(allData, "./ailMasterData.txt", sep="\t", row.names= FALSE, quote=FALSE)
covariates <- allData[c(1:36, 110:115, 120:123, 131:132, 134:135, 138:141)]
covariates <- covariates[-c(45)]
names(covariates)
phenotypes2 <- allData[-c(2:24, 26,28,30,32,34,36,  
                         110:115, 121:123, 131, 134, 138)]


write.table(covariates, "./ailCovariates.txt", sep="\t", row.names= FALSE, quote=FALSE)
write.table(phenotypes2, "./ailPhenotypes.txt", sep="\t", row.names= FALSE, quote=FALSE)

# went in and switched col headers around
phen <- read.table("./ail.txt", sep="\t", header=T)
cov <- read.table("./ailCovariates.txt", sep="\t", header=T)

phenoColNames <- names(phen)
covColNames <- names(cov)

write.table(covColNames, "./covColNames.txt", sep="\t", quote=FALSE)
write.table(phenoColNames, "./phenColNames.txt", sep="\t", quote=FALSE)


## info for Ari
ariData <- allData[c(1:4, 7,8,13,14,18,22,134,135)]
names(ariData)
head(ariData)
write.table(ariData, "./AIL-hindlimb-sampleInfo.txt", sep="\t", row.names=FALSE,quote=FALSE)
