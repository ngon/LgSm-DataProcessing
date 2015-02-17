setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/pedigree")
ped=read.table("ailpedigree.csv", sep=",", header=TRUE)

# What line is the first mouse from F45 in?
ped[ped$id == "39369.1",]
ped[6700:6715,]

# Trim ped to include only F45-56
pedtrim <- ped[6706:10173,]
# remember to add --no-pheno when running ped

write.table(file="ailped.45to56.ped", sep="\t", 
            row.names=FALSE, col.names=FALSE, pedtrim)


__________________________________

phenos <- read.table("allData.txt", sep="\t", header=T)
phenos[phenos$dam == "50267",]

___________________________________

## preparing files for pedchecking in Mark Abney's idcf software

setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/pedigree")

## create pedigree file ##

ped=read.table("ailpedigree.csv", sep=",", header=TRUE, as.is=TRUE)
f34ped <- read.table("pedigreeF34.csv", sep=" ", header=TRUE, as.is=TRUE)

f2s <- f34ped[f34ped$generation == "F2",]

# these F2 animals are listed in the AI_Line info F3 file rather than in the F2 file.
# therefore they are out of order when the ped is generated in R.
# i moved them to the correct lines (232-4) so I could replace their sire and dam ids with the corresponding ids in the F34 pedigree file (i.e. the pedigree file that was used in the 2009 Genetics paper). This file codes F1 parents as -9997 and -9996, and includes the F1 and F0 generations. 

# mice that were out of order:
ped[ped$id == "1182",]
ped[ped$id == "1190",]
ped[ped$id == "1191",]

f34ped <- cbind(f34ped$id, f34ped$sire, f34ped$dam)
ped <- cbind(ped$id, ped$sire, ped$dam)
ped <- ped[-c(1:234),]
ped <- rbind(f34ped[1:238,], ped)
head(ped)

# change F0 founders to id "0"
ped[1,1] <- 0
ped[2,1] <- 0 
ped[3, 2:3] <- 0
ped[4, 2:3] <- 0
head(ped)
write.table(ped, file="ailpedforidcf.txt", sep=" ", row.names=FALSE, col.names=FALSE)


## create study file ##

setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles")
testmice <- read.table(file="allData.txt", sep="\t", header=TRUE)

study <- as.data.frame(testmice$id)
study <- study[1] - 0.1
head(study)
write.table(study, file="studyfileidcf.txt", sep=" ", row.names=FALSE, col.names=FALSE)

________________

# idcf output spit out the following. These are animals from Cheverud generations. 

# Mother of 41095 not found in pedigree.
# Father of 41095 not found in pedigree.
# Father of 50924 not found in pedigree.
# Mother of 51914 not found in pedigree.
# Father of 51914 not found in pedigree.
# Mother of 52016 not found in pedigree.
# Father of 52016 not found in pedigree.


pedigree=read.table("ailpedigree.csv", sep=",", header=TRUE, as.is=TRUE)

# changing these unlisted parents to "0" for this purpose

ped2 <- as.data.frame(ped)
names(ped2)[1:3] <- c("id", "sire", "dam")

ped2[ped2$id == "41095",]
ped2[3212, 2:3]<- 0

ped2[ped2$id == "50924",]
ped2[3867, 2:3]<- 0

ped2[ped2$id == "51914",]
ped2[4004, 2:3]<- 0

ped2[ped2$id == "52016",]
ped2[4005, 2:3]<- 0


write.table(ped, file="ailpedforidcf.txt", sep=" ", row.names=FALSE, col.names=FALSE)


### update: idcoef doesn't accept ids unless they're integers. i replaced all instances of ".1" with "0" in the ped file. 0.1 was added to ids from F33 and beyond to avoid overlap with ids from Cheverud generations. all Cheverud ids are 5 digits or less, so I added additional zeros to make all ids form Palmer generations 6 digits long. Some are 7 digits because they already had an additonal zero added (e.g. to correct errors - see ped readme file). There were also two IDs from Cheverud generations that were inexplicably labeled "50.219" and "51.224". I changed these to "5021900" and "5122400", respectively. 

# testing to make sure it looks ok after saving in excel text format
check <- read.table(file="ailpedforidcf.txt", sep="\t")
head(check)

# remove extra columns
check <- check[1:3]
head(check)
write.table(check, file="ailped.integer.txt", sep=" ", row.names=F, col.names=F)






