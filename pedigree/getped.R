setwd("C:/Users/Administrator/Desktop/AILped")
source("kinship.r")
source("find_mates.r")
data.dir = "C:/Users/Administrator/Desktop/AILped/"


### get ped for AIL breeders ###

ped = c()

for (i in 2:56){
  file.name = paste(data.dir,"AI Line Info_",i,".csv",sep="")
  ped.tmp = read.table(file=file.name, sep=",", skip = 1,
     as.is=TRUE, na.strings="?")
  ped.tmp = ped.tmp[,c(1,7,8,3)] ## this differs from original breeder code in that
                                 ## cols 7 and 8 are switched to match PLINK format
  dim(ped.tmp)
  if (i == 2){
    ped.tmp[,c(2,3)] = 0}
 
 if (i == 34){
    ## Determine which ids already have decimal; add 0.1 to the rest##
    idx = (ped.tmp[,1] - ped.tmp[,1]%/%1) == 0
    ped.tmp[idx,1] = ped.tmp[idx,1]+.1 }

  # CHECK THAT PARENTS ARE IN PREVIOUS GENERATION #
  # see Readme file in data.dir for details on output # 
  if (i != 2){
    tmp = check.ped(ped.tmp , ped.prev)
    if (length(tmp) > 0){
       print(paste("Generation:",i,"Parents not present in previous generation:",tmp))
       }
   }
    ped.prev = ped.tmp	
  ped = rbind(ped,ped.tmp)
}

# format sex coding so that M=1 and F=2 
ped[,4] = as.character(ped[,4])
ped[ped[,4]=="F",4] = 2
ped[ped[,4]=="M",4] = 1
ped[,4] = as.numeric(ped[,4])

# name cols
colnames(ped) <- c("id", "sire", "dam", "sex")
fam <- rep(1, length(ped$id))
ped <- cbind(fam, ped)

# # check for duplicates
# brdups <- duplicated(ped$id)
# sum(brdups)


### reformat ped info for F50-56 AILs that were tested/genotyped ###

# get data
testmice<- read.table(file="./allData.txt", sep="\t", header=T)
testped <- cbind(testmice$id, testmice$sire, testmice$dam, testmice$sex)
testped <- as.data.frame(testped)

# change sex coding from M=0 and F=1 to M=1 and F=2
testped$V4 <- as.character(testped$V4)
testped[testped[,4]=="2",4] = 0
testped[testped[,4]=="1",4] = 2
testped[testped[,4]=="0",4] = 1

# name cols
colnames(testped) <- c("id", "sire", "dam", "sex")
fam <- rep(1, length(testped$id))
testped <- cbind(fam, testped)

# reformat ids to match the breeder ped
testped$id <- testped$id + 0.1
testped$sire <- testped$sire + 0.1
testped$dam <- testped$dam + 0.1

# look for duplicates
testdups <- duplicated(testped$id)
sum(testdups)

# combine ail and breeder peds
fullped <- rbind(ped, testped)
write.table(file="./ailpedigree.csv", sep=",", row.names=FALSE, fullped)
write.table(file="./ailped.F49to56.csv", sep=",", row.names=FALSE, fullped[7435:10173,])

# look for duplicate ids
# dups <- duplicated(fullped$id)
# dups <- which(dups==1)
# fullped[dups,]


#### Fix my phenotype files so that they have the corrected dam and sire ids

setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles")
cov <- read.table(file="./covariates-edit dam sire.txt", sep="\t", header=T)
covariateCols <- (names(cov))
cpp <- read.table(file="./cppData.all-edit dam sire.txt", sep="\t", header=T)
methcols <- names(cpp)

data <- read.table(file="./allData.txt", sep="\t", header=T)
dataCols <- names(data)

covCols=which(dataCols %in% covariateCols)
covariates=data[,covCols]
covariates = cbind(cov[,1], covariates, cov[,46])
names(covariates)[1] <- "one"
names(covariates)[46] <- "coat"
write.table(file="./covariates.txt", sep="\t", covariates)

cppCols <- which(dataCols %in% methcols)
cppPheno <- data[,cppCols]
write.table(file="./cppPhenotypes.txt", sep="\t", cppPheno)

