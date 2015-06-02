### PURPOSE: Check for mismatches in phenotype and genotype data for F50-56

### LOAD DATA ------------------------------------------------------------------
### geno.samples is genotyped samples, info is phenotype data with fam info
### ail.kinship is a matrix of kinship coeffs from the pedigree,and grm is the
### IBS/grm for all genotyped samples.
geno.samples <- read.table("./genotyped.samples.txt",as.is=T)
names(geno.samples) <- c("id")
geno.samples$id <- as.numeric(geno.samples$id)
geno.samples$id <- geno.samples$id + 0.1

# all phenotype and covariate data
info <- read.table("./dataFiles/phenotypes/allData.txt", sep="\t",
                   header=T, as.is=T)[1:6]
info$id <- info$id + 0.1
info$cc=NULL

#which(!info$id %in% row.names(ail.kinship)) # 89, 46368, gen 50
#info <- info[c(1:88, 90:1124),]


# kinship and grm matrices
load("./pedigree/ail.kinship.QTLRel.RData")
grm <- read.table("./dataFiles/chrAll.cXX.txt", as.is=TRUE)
row.names(grm)[1:1830] <- geno.samples$id[1:1830]
names(grm)[1:1830] <- geno.samples$id[1:1830]

### SIBLING PAIRS & PED DATA ---------------------------------------------------
### Is ped kinship between sibs (defined as having the same fam name) around
### the expected value of 0.5
### ANSWER: Min kinship for all sibs is 0.4128; max is 0.5763.
### There are two sib pairs with k < 0.5:
### 46023 & 46024, different parents (but both listed as fam BrF49-24).
      id gen      fam   dam  sire
46023.1  50 BrF49-24 45383 45318
46024.1  50 BrF49-24 45392 45446 # should be BrF49-28
### 50210 & 50377. same as above, both listed as BrF51-34
      id gen      fam   dam  sire
50210.1  52 BrF51-34 49017 48404 # is a sib with 50377. change parents 48442 and 49096, as below.
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
### we have no misnatches.
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

sireNames <- info$sire[which(duplicated(info$sire))]
sibAa <- c()
sibBb <- c()
kinship.s <- c()
for (sire in sireNames){
    sibs <- info[info$sire == sire,]
    sibs$id <- as.character(sibs$id)
    sibAa <- append(sibAa, sibs[1,1])
    sibBb <- append(sibBb, sibs[2,1])
    kinship.s <- append(kinship.s, ail.kinship[sibs[1,1], sibs[2,1]])
}
siblings.s <- data.frame(sibAa, sibBb, kinship.s)

### filtering out sibs by sire or dam gives different results. there are 542
### dams that are parents of siblings, yet 546 sires that are parents of siblings.
dim(siblings.d); dim(siblings.s)
### there are 581 unique dams, yet only 577 unique sires. the nums should be equal.

stest<- info[which(info$sire %in% sireNames),]
dtest<- info[which(info$dam %in% damNames),]

info[which(!sibA %in% sibAa),]

         id gen      fam   dam  sire
127 46490.1  50 BrF49-74 45347 45470 # ok
128 46491.1  50 BrF49-74 45347 45470 # ok
131 46640.1  50 BrF49-70 45437 45470 # dam, sire = 45424, 45394
132 46641.1  50 BrF49-70 45437 45470 # "
         id gen      fam   dam  sire
390 50225.1  52 BrF51-56 49049 49500 # sire = 49000
392 50230.1  52 BrF51-56 49049 49500 # sire = 49000
491 50716.1  53 BrF52-28 50219 49500 # this should be 50714. 50716 is his sib, a breeder mouse. since they're brothers it really doesn't matter and i'm leaving it as is for now.it looks like their ID tags were switched.
494 50723.1  53 BrF52-28 50219 49500 # 49500 is an F52 mouse.
         id gen      fam   dam  sire
536 51211.1  53 BrF52-60 50393 50138 # sire = 50124
539 51220.1  53 BrF52-60 50393 50138 # sire = 50124
572 51289.1  53 BrF52-61 50389 50138 # ok
574 51294.1  53 BrF52-61 50389 50138 # ok
          id gen      fam   dam  sire
981  54282.1  55 BrF54-77 51957 52184 # ok
986  54291.1  55 BrF54-77 51957 52184 # ok
1060 54853.1  56 BrF55-12 52097 52184 # sire = 52084
1061 54854.1  56 BrF55-12 52097 52184 # sire = 52084



### SIBLING PAIRS AND GRM ------------------------------------------------------
### find families with more than 1 sibling (should only be a few)
levs <- as.factor(info$fam)
familySizes <- which(table(levs) > 2)
BrF50-17 BrF50-82 BrF51-10 BrF52-57 BrF52-65 BrF54-51 BrF55-04
(3, 4, 3, 3, 3, 3, 4)

# look for typos
all.not <- geno.samples[!geno.samples$id %in% alldata$id,]
gen<-c(42,42,42,42, 43, 50,50,50,50,50,50,50,50,50,50,50,50,50,51,51,51,51,51,NA,51, 51, 53,53,52,52, 52,52,53,53,53,53,53,54,54,54)
sex<- c(1,1,1,1, 2, 1, 2,2,2,1,1,1,1,1,1,2,2,1,2,2,1,1,2, NA, 2,1, 1,2, 2,1,2,2,2,2,1,1,2,1,1,2)
missingRows <- cbind(as.integer(all.not[1:40]), as.integer(gen), as.integer(sex))



which(grm['49290.1',]>0.05) # sib 49289. gen 51. AIL lib 42. this one looks like a typo in the testing/breeding spreadsheet. it's entered as 48290 even though all the surrounding ids are in the 49000s. so althought 49290 is probably the correct id, i'm leaving it as 48290.

"54386" # flowcell 28 lib 84 - gen 56
which(grm['54386.1',] > 0.05) # 54385 is definitely a sib. 52197 could be a sib, or is at least related. gr between 85,86 = 0.097; with 52197 = 0.051

"57801" # flowcell 25 lib 65 - allegedly generation 54 - should be 51801.
which(grm['57801.1',] > 0.05)
# kinship from the GRM adds up - if 51801, i'd expect her to be most related to the potential sibi=ling, 51794, and that is the case.

relatives<- c()
for (i in rownames(grm)){
    #if(length(which(grm[i,] > 0.05)) > 1){
        relatives[[i]] <- which(grm[i,] > 0.05)
       # }
}

r <- names(relatives)
r <- as.integer(r)   # list of 1076


## kinship coeffs
library(QTLRel)
ail.ped <- read.table("./pedigree/pedforQTLRel.txt", sep="\t", header=T)[1:5]
ail.k <- kinship(ped=ail.ped, ids=geno.samples$id)
save(ail.k, file="./pedigree/kinship.geno.samples.RData")


check.kinship <- function(ailrel){
values <- c()
for (name in seq_along(names(ailrel))){
    len <- length(ailrel[[name]])
    id.tmp <- as.numeric(names(ailrel)[name]) - 0.1
    id <- ailrel[[name]][[which(names(ailrel[[name]]) == id.tmp)]]
    values[[name]] <- data.frame(idx1=rep(id, times=len),
                                 idx2=ailrel[[name]])
    k.grm <- c()
    k.ped <- c()
    ped.kin <- c()
        for (i in seq_along(values[[name]][['idx2']])){
           k.grm <- c(k.grm, grm[values[[name]][1,1], values[[name]][i,2] ])
           k.ped <- c(k.ped, ail.k[values[[name]][1,1], values[[name]][i,2] ])
           ped.kin <-c(ped.kin, paste(names(which(ail.k[values[[name]][i,2],] > 0.48)),
                            sep="", collapse=" "))
        }
    rownames(values[[name]]) <- as.integer(rownames(values[[name]]))+0.1

    values[[name]] <- cbind(values[[name]],
                            ail.ped[ail.ped$id %in% rownames(values[[name]]), c(2:5)])

    values[[name]][["k.grm"]] <- k.grm
    values[[name]][["k.ped"]] <- k.ped
    values[[name]][["ped.kin"]] <- ped.kin
}
names(values)<- names(ailrel)
return(values)
}

## scaled and centered values
tk <- scale(ail.k)
t1 <- scale(grm, center=FALSE, scale=apply(grm, 2, sd, na.rm=T))
ail.k <- tk
grm <- t1


kinCheck.allgenoScaled<- check.kinship(ailrel=relatives)
save(kinCheck.allgenoScaled, file="kinCheck.allgenoScaled.RData")
#save(ail.k,file="./pedigree/kinship.geno.samples.RData" )
#save(grm, file="./dataFiles/grm.geno.samples.RData")

df<-c()
for(i in kinCheck.allgenoScaled){
    df <- rbind.data.frame(df,i)
}

df$dam <- as.factor(df$dam)
df$generation <- as.factor(df$generation)


## unscaled/uncentered values
# sibMean.grm<- tapply(df$k.grm, df$dam, mean)
# sibMean.ped<- tapply(df$k.ped, df$dam, mean)
# sibMed.grm<- tapply(df$k.grm, df$dam, median)
# sibMed.ped<- tapply(df$k.ped, df$dam, median)
# sibStats <- cbind(mean(sibMean.grm), mean(sibMean.ped), mean(sibMed.grm),
#                   mean(sibMed.ped), mean(grm), mean(ail.k),
#                   median(grm), median(ail.k))
> grm<- grm_scaled + 1
> ail.k <- ail.k_scaled +1

grmk <- rowMeans(grm)
pedk <- rowMeans(ail.k)
idk <- rownames(grm)
relAll <- cbind(idk, grmk, pedk, ail.ped[ail.ped$id %in% idk,c(4:5)])
gMeanAll<- tapply(relAll$grmk, relAll$generation, mean)
pMeanAll<- tapply(relAll$pedk, relAll$generation, mean)
allStats <- rbind(gMeanAll, pMeanAll)
allStats <- data.frame(allStats[,c(34, 36:45, 47:53)])
allStats$stat <- c("grm.mean", "ped.mean")
allStats <- melt(allStats, id.vars="stat")
names(allStats)[2] <- "generation"

sans.pedk <- rowMeans(pedsansRel)
sans.grmk <- rowMeans(sansRel)
idsansk <- rownames(sansRel)
relSansSibs <- cbind(idsansk, sans.grmk, sans.pedk, ail.ped[ail.ped$id %in% idsansk, c(4:5)])
gMeanSans<- tapply(relSansSibs$sans.grmk, relSansSibs$generation, mean)
pMeanSans<- tapply(relSansSibs$sans.pedk, relSansSibs$generation, mean)
sansStats <- rbind(gMeanSans, pMeanSans)
sansStats <- data.frame(sansStats[,c(34, 36:45, 47:53)])
sansStats$stat <- c("grm.mean", "ped.mean")
sansStats <- melt(sansStats, id.vars="stat")
names(sansStats)[2] <- "generation"
# cleanup
rm(gMeanAll, pMeanAll, grmk, pedk, sans.pedk, sans.grmk, idk, idsansk)



p<-
    ggplot(data=allStats, aes(x=generation))+
    geom_bar(data=subset(allStats, stat == 'grm.mean'),
             aes(y=value),
             stat='identity') +
    geom_bar(data=subset(allStats, stat == 'ped.mean'),
             aes(y= -value),
             stat='identity') +
    scale_y_continuous(limits=c(-1.4, 1.2),
                       breaks=c(seq(-1.4, 1.2, 0.1)))+


                       #breaks=c(-35, -25, -24.3, -15,-14.94, -5, 5, 14.66, 15, 25, 30,
                                #31.2, 35),
                       #labels=c(35, " ", paste0('self ', 24.3), " ",
#                                 paste0('sib ',14.9), 5, 5, paste0('sib ',14.6)," ",
#                                 25, " ", paste0('self ',31.2),35))+

    xlab("Generation")+
    ylab("Kinship\n Pedigree - Genotypes")+
    theme_bw() +

    geom_hline(yintercept=0, color='white')


    geom_hline(yintercept=31.2, color='dodgerblue', lty=2)+
    geom_hline(yintercept=-24.3, color='dodgerblue', lty=2)+
    geom_hline(yintercept=c(14.66, -14.94), color='black', lty=2)+

    scale_fill_gradient2(low="dodgerblue", high="dodgerblue", mid="white",
                         midpoint=0) +
    theme(legend.position="none")







write.table(df, "./dataFiles/kinCheck.allgeno.txt", sep=" ", col.names=T,
            row.names=F, quote=F)

p+
geom_bar(data=subset(sansStats, stat == 'grm.mean'),
         aes(y=value, alpha=0.5),
         stat='identity') +
    geom_bar(data=subset(sansStats, stat == 'ped.mean'),
             aes(y= -value, alpha=0.5),
             stat='identity')

geom_hline(yintercept=0, color='white')+
scale_fill_gradient2(low="dodgerblue", high="dodgerblue", mid="white",
                     midpoint=0) +
    theme(legend.position="none")


### FUNCTION: ALL DUP
### allDup finds and keeps all duplicated records.
allDup <- function (x){duplicated(x) | duplicated(x, fromLast = TRUE)}
### list of mice with no siblings
#sibMice <- ail.ped[allDup(ail.ped$dam),] # 1070
#soloMice <- ail.ped[!allDup(ail.ped$dam),] # 53
sansSibs <- ail.ped[which(duplicated(ail.ped$dam) == FALSE),] # 581

#sibRows <- which(rownames(grm) %in% sibMice$id)
#soloRows <- which(rownames(grm) %in% soloMice$id)
sansRows <- which(rownames(grm) %in% sansSibs$id)

#sibRel <- grm[sibRows, sibRows]
#soloRel <- grm[soloRows,soloRows]
sansRel <- grm[sansRows, sansRows]
pedsansRel <- ail.k[sansRows,sansRows]

## unscaled/uncentered values

genMean.grm<- tapply(df$k.grm, df$generation, mean)
genMean.ped<- tapply(df$k.ped, df$generation, mean)
genMed.grm<- tapply(df$k.grm, df$generation, median)
genMed.ped<- tapply(df$k.ped, df$generation, median)

genStats <- rbind(genMean.grm, genMean.ped, genMed.grm, genMed.ped)
genStats <- data.frame(genStats[,c(34, 36:45, 47:53)])
genStats$type <- c("grm.mean", "ped.mean", "grm.median", "ped.median")


# library(lattice)
# test<-colorRampPalette(cols, alpha = TRUE)(50)
# cols<- brewer.pal(11, 'PiYG')
# levelplot(soloRel, col.regions=test, cuts=49, cex.axis=0.8,
#           scales=list(x=list(rot=45)))
# levelplot(sibRel[1:20,1:20], col.regions=test, cuts=49, cex.axis=0.8,
#           scales=list(x=list(rot=45)))
# levelplot(sansRel[1:20,1:20], col.regions=test, cuts=49, cex.axis=0.8,
#           scales=list(x=list(rot=45)))

