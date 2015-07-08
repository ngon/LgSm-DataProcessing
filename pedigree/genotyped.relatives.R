### PURPOSE: MAKE A LIST OF SIB PAIRS GENOTYPED FROM F50-55, AND MAKE A LIST OF
### TRIOS/FAMILIES GENOTYPED FROM F39-43. THESE CAN BE USED FOR ERROR CHECKING
### FURTHER DOWN THE LINE. 

### THIS WAS USED TO CREATE /PEDIGREE/GENOTYPEDRELATIVES.R


load("./pedigree/info.phenoSamples.RData")
info <- info[order(info$gen, info$dam),]
# get info for genotyped samples only (ids)
info$id <- info$id - 0.1
info <- info[info$id %in% ids,]

sibPairs <- tapply(info$id, INDEX=info$dam, FUN=as.vector, simplify=FALSE)
for (i in seq_along(sibPairs)) {
    if (length(sibPairs[[i]]) < 2) {
        sibPairs[i] <- NULL
    }
}
# MORE VERBOSE THAN TAPPLY
# retain only mice with siblings
# sibs.tmp <- c()
# sibs <- list()
# for (i in 1:(length(info$dam)-1)) {
#     if (info$dam[i] == info$dam[i+1]) {
#         sibs.tmp <- rbind(sibs.tmp, info[i,], info[i+1,])
#     } # gets df containing sib pair info
#     if (info$dam[i] == info$dam[i+1]) {
#         sibs[[i]] <- c(info[i,1], info[i+1,1])
#     } # makes list of sibs for each dam
# }
# sibs <- sibs[!sapply(sibs, is.null)]
# names(sibs) <- unique(sibs.tmp$dam)


# get trios for F39-43
ped <- read.table("./pedigree/pedforQTLRel.txt", sep="\t", header=T)[1:5]
ids <- ids + 0.1
ped <- ped[ped$id %in% ids,]
sires <- ped[ped$sire %in% ids,] 
dams <- ped[ped$dam %in% ids,] 

trios <- dams[dams$id %in% sires$id,]
trios[c(1,2,3)] <- trios[c(1,2,3)] - 0.1
trios <- trios[1:3]
trios <- trios[order(trios$sire),]
fams <- tapply(trios$id, INDEX=trios$fams, FUN=as.vector, simplify=FALSE)
sires2 <- tapply(trios$sire, INDEX=trios$sire, FUN=as.vector, simplify=FALSE)
dams2 <- tapply(trios$dam, INDEX=trios$dam, FUN=as.vector, simplify=FALSE)

# sires to append
siresOrdered <- c()
for (i in seq_along(sires2)) {
    siresOrdered <- c(siresOrdered, sires2[[i]][1])
}
# dams to append
damsOrdered <- c()
for (i in seq_along(dams2)) {
   damsOrdered <- c(damsOrdered, dams2[[i]][1])
}
# list of families with all genotyped members
for (i in seq_along(fams)) {
    fams[[i]] <- append(fams[[i]], c(siresOrdered[i], damsOrdered[i])) 
}

# save(fams, sibPairs, file="./pedigree/genotypedRelatives.RData")
