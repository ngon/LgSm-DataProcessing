### PURPOSE: Check for mismatches in phenotype and genotype data.

### LOAD DATA ------------------------------------------------------------------
### geno.samples is genotyped samples, info is phenotype data with fam info
### ail.kinship is a matrix of kinship coeffs from the pedigree,and grm is the
### IBS/grm for all genotyped samples.
geno.samples <- read.table("./genotyped.samples.txt",as.is=T)
geno.samples[which(geno.samples$V1 == "54109_11"), 1] <- "54109"
geno.samples[which(geno.samples$V1 == "51717-19"), 1] <- "51717"
geno.samples[which(geno.samples$V1 == "51348"), 1]    <- "51349"
names(geno.samples) <- c("id")
geno.samples$id <- as.numeric(geno.samples$id)
geno.samples$id <- geno.samples$id + 0.1

info <- read.table("./dataFiles/phenotypes/allData.txt", sep="\t",
                   header=T, as.is=T)[1:6]
info$id <- info$id + 0.1
info$cc=NULL

load("./pedigree/ail.kinship.QTLRel.RData")
grm <- read.table("./dataFiles/chrAll.cXX.txt", as.is=TRUE)
row.names(grm)[1:1830] <- geno.samples$id[1:1830]
names(grm)[1:1830] <- geno.samples$id[1:1830]

### SIBLING PAIRS & PED DATA ---------------------------------------------------
### Is ped kinship between sibs (defined as having the same fam name) around
### the expected value of 0.5?
### ANSWER: Min kinship for all sibs is 0.4128; max is 0.5763.
### There are two sib pairs with k < 0.5:
### 46023 & 46024, different parents (but both listed as fam BrF49-24).
id gen      fam   dam  sire
46023.1  50 BrF49-24 45383 45318
46024.1  50 BrF49-24 45392 45446 # should be BrF49-28
### 50210 & 50377. same as above, both listed as BrF51-34
id gen      fam   dam  sire
50210.1  52 BrF51-34 49017 48404 # should be BrF51-47
50377.1  52 BrF51-34 48442 49096
### The dam and sire ids are CORRECT in both cases; only the fam name is wrong.

# There are 582 unique families (541 sib pairs = 1082 mice; 41 w/o sib)
famNames <- info$fam[which(duplicated(info$fam))]
sib1 <- c()
sib2 <- c()
kinship <- c()
for (family in famNames){
    sibs <- info[info$fam == family,]
    sibs$id <- as.character(sibs$id)
    sib1 <- append(sib1, sibs[1,1])
    sib2 <- append(sib2, sibs[2,1])
    kinship <- append(kinship, ail.kinship[sibs[1,1], sibs[2,1]])
    }
siblings <- data.frame(sib1, sib2, kinship)

### Same as above, except sibs are paired by unique dam names instead. Here
### we have no misnatches. Perhaps throw out the fam.id altogether; use parent
### IDs as substitutes for fam.id.
damNames <- info$dam[which(duplicated(info$dam))]
sibA <- c()
sibB <- c()
kinship.d <- c()
for (dam in damNames){
    sibs <- info[info$dam == dam,]
    sibs$id <- as.character(sibs$id)
    sibA <- append(sibA, sibs[1,1])
    sibB <- append(sibB, sibs[2,1])
    kinship.d <- append(kinship.d, ail.kinship[sibs[1,1], sibs[2,1]])
}
siblings.d <- data.frame(sibA, sibB, kinship.d)

### SIBLING PAIRS AND GRM ------------------------------------------------------

### Make a new info df for AILs that were phenotyped and genotyped (29 mice
### remain to be genotyped)
grinfo <- info[info$id %in% row.names(grm),]


sampDams<- grinfo$dam[which(duplicated(grinfo$dam))]
sib1 <- c()
sib2 <- c()
sib3 <- c()
ibs <- c()
for (dam in sampDams){
    sibs <- grinfo[grinfo$dam == dam,]
    sibs$id <- as.character(sibs$id)
    sib1 <- append(sib1, sibs[1,1])
    if(length(sibs$id == 3)){
        sib2 <- append(sib2, sibs[2,1])
        sib3 <- append(sib3, sibs[3,1])
    } else if(length(sibs$id == 2)){
        sib2 <- append(sib2, sibs[2,1])
        sib3 <- append(sib3, "NA")
    } else if(length(sibs$id == 1)){
        sib2 <- append(sib2, "NA")
        sib3 <- append(sib3, "NA")
    }

    ibs <- append(ibs, grm[sibs[1,1], sibs[2,1]])
}
sibs.ibs <- data.frame(sib1, sib2, ibs, sib3)

### find families with more than 1 sibling (should only be a few)
levs <- as.factor(info$fam)
familySizes <- which(table(levs) > 2)
BrF50-17 BrF50-82 BrF51-10 BrF52-57 BrF52-65 BrF54-51 BrF55-04
(3, 4, 3, 3, 3, 3, 4)

### FUNCTION: ALL DUP
### allDup finds and keeps all duplicated records.
allDup <- function (x){duplicated(x) | duplicated(x, fromLast = TRUE)}
### list of mice with no siblings
sibMice <- info[allDup(info$dam),] # 1070
soloMice <- info[!allDup(info$dam),] # 53
sansSibs <- info[which(duplicated(info$dam) == FALSE),] # 581

sibRows <- which(rownames(grm) %in% sibMice$id) # 1042
soloRows <- which(rownames(grm) %in% soloMice$id) # 52
sansRows <- which(rownames(grm) %in% sansSibs$id) # 565

sibRel <- grm[sibRows, sibRows]
soloRel <- grm[soloRows,soloRows]
sansRel <- grm[sansRows, sansRows]

library(lattice)

test<-colorRampPalette(cols, alpha = TRUE)(50)
cols<- brewer.pal(11, 'PiYG')


levelplot(soloRel, col.regions=test, cuts=49, cex.axis=0.8,
          scales=list(x=list(rot=45)))

levelplot(sibRel[1:20,1:20], col.regions=test, cuts=49, cex.axis=0.8,
          scales=list(x=list(rot=45)))

levelplot(sansRel[1:20,1:20], col.regions=test, cuts=49, cex.axis=0.8,
          scales=list(x=list(rot=45)))

