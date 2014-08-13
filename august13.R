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

allData <- merge(cppData, ppiOther, by="id", all=TRUE)
write.table(allData, "./mergedData.txt", sep="\t", row.names=FALSE, quote=FALSE)

# sex, fam/sire/dam, and cage variables x and y do not match in all cases
# this is because in one dfs the variable was missing, and in the other it was not
check <- allData[complete.cases(allData),]
