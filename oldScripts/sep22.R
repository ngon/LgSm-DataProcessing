## Format phenotype files (complete)

list.files()

### Read data

######## TO RUN WHEN DISSECTION IS COMPLETE ##############
# physiological phenos 
# other <- read.table("./otherPhenos.ALL.txt", sep="\t", na.strings=NA, header=TRUE)
# other <- other[do.call(order, other),]
########################################################

# cpp data
cppData<- read.table("./dataFiles/cppData.all.txt", 
                     sep="\t", header=TRUE, na.strings="NA")
# read ppi info 
ppiInfo <- read.table(file="./dataFiles/ppi.info.all.txt", 
                      sep="\t", header=TRUE, na.strings="NA")
# read ppi data
ppiData <- read.table(file="./dataFiles/ppi.data.ALL.txt", 
                      sep="\t", header=TRUE, na.strings="NA")
# read other phenos
other <- read.table(file="./dataFiles/otherPhenos.all.txt",
                        sep="\t", header=TRUE, na.strings="NA")


# sort data frames
cppData <- cppData[do.call(order,cppData),]
ppiInfo <- ppiInfo[do.call(order, ppiInfo),]
ppiData <- ppiData[do.call(order, ppiData),]
other <- other[do.call(order, other),]

#### format ppi data
# pull out relevant columns from ppiData and rename some of them
ppiData <- ppiData[,c("ID", "BOX","PPI3", "PPI6", "PPI12",
                      "avg_ppi", "startle", "habituation",
                      "sum.nostim")]
names(ppiData)[1:6] <- c("id", "ppi.box", "ppi3", 
                         "ppi6", "ppi12", "avg.ppi")

# merge ppiData and ppiInfo
ppiMerge <- merge(ppiInfo, ppiData, by="id", all=TRUE)

#check for and remove duplicate values (there is 1)
which(duplicated(ppiMerge$id)) 
ppiMerge[1002:1005,] 
ppiMerge <-ppiMerge[-c(1004),]# 53420 - remove entry with ppi.weight = NA

# remove blank rows in otherData
other <- other[-c(1:3),]

# merge cpp and other
cppMerge <- merge(cppData, other, by="id", all=TRUE)
which(duplicated(cppMerge$id))
cppMerge[807:808,] 
cppMerge<-cppMerge[-c(808),] # remove incorrect id

# wrote file and opened it. removed a few cells with typos and saved the file as cppCheck.
# cppCheck.txt is updated and can be loaded without running the steps above. 
write.table(cppMerge, "./dataFiles/cppCheck.txt", sep="\t", row.names=F, quote=F)

cppMerge<- read.table("./dataFiles/cppCheck.txt", sep="\t", header=T, na.strings="NA")
# remove duplicate cols because they are identical
cppMerge <- cppMerge[c(1:127)]


# merge all pheno data
ppiMerge$id <- as.factor(ppiMerge$id)
allData <- merge(cppMerge, ppiMerge, by="id", all=TRUE)
which(duplicated(allData$id)) # no duplicates. good. 
allData <- allData[do.call(order, allData),]

# remove rows with NAs (pheno errors from early generations; IDs with hyphens or slashes)
# these mice will have NA under allData$sex; all other mice have either M or F. 
blankValues <- which(is.na(allData$sex))
# remove rows with allData$sex == NA
allData <- allData[-c(blankValues),]

# rename a few columns
names(allData)[120:127]<- names(other[2:9])

# make separate coat color categories
unique(allData$cc)
which(allData$cc == "AGR/L") # 587
allData[587, 3] <- "AGRL"
which(allData$cc == "ARG") #535
allData[535, 3] <- "AGR"
# fill in some missing CC values
allData[95,1:3]
allData[95,3] <- "BL"
allData[93,1:3]
allData[93,3] <- "AGR"
allData$cc <- droplevels(allData$cc)

# colors were inconsistently entered. good to have the ear tag info, but
# better to have a separate variable for color alone. 
colors <- allData$cc
levels(colors)[levels(colors)=="AGL"] <- "A"
levels(colors)[levels(colors)=="AGR"] <- "A"
levels(colors)[levels(colors)=="AGRL"] <- "A"
levels(colors)[levels(colors)=="AL"] <- "A"
levels(colors)[levels(colors)=="AR"] <- "A"
levels(colors)[levels(colors)=="ARL"] <- "A"
levels(colors)[levels(colors)=="WL"] <- "W"
levels(colors)[levels(colors)=="WR"] <- "W"
levels(colors)[levels(colors)=="WRL"] <- "W"
levels(colors)[levels(colors)=="BL"] <- "B"
levels(colors)[levels(colors)=="BR"] <- "B"
levels(colors)[levels(colors)=="BRL"] <- "B"


# save all pheno data in one master file
write.table(allData, "./dataFiles/allData.txt", sep="\t", row.names=FALSE, quote=FALSE)

covariates.all <- allData[c(1:2, 4:36,120,121,123,124,126:130)]
one <- rep(1, 1123)
covariates.all<- cbind(one, covariates.all, colors)
names(covariates.all)[46] <-"coat"

phenos.all <- allData[c(1, 109,119, 44:51, 94:101, 52:58, 102:108,
                       37:43, 59:93,
                       131:136, 122, 125, 126)]


write.table(covariates.all, "./dataFiles/covariates.txt", sep="\t", row.names= FALSE, quote=FALSE)
write.table(phenos.all, "./dataFiles/phenotypes.txt", sep="\t", row.names= FALSE, quote=FALSE)

write.table(covariates.all, "./dataFiles/cov.noHeader.txt", sep="\t", row.names= FALSE, col.names=FALSE, quote=FALSE)
write.table(phenos.all, "./dataFiles/pheno.noHeader.txt", sep="\t", row.names= FALSE, quote=FALSE, col.names=FALSE)


phenoColNames <- names(phenos.all)
covColNames <- names(covariates.all)

write.table(covColNames, "./covColNames.txt", sep="\t", quote=FALSE)
write.table(phenoColNames, "./phenColNames.txt", sep="\t", quote=FALSE)

# 
# ## info for Ari
# ariData <- allData[c(1:4, 7,8,13,14,18,22,134,135)]
# names(ariData)
# head(ariData)
# write.table(ariData, "./AIL-hindlimb-sampleInfo.txt", sep="\t", row.names=FALSE,quote=FALSE)
