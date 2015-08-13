load("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/pedigree/grm.geno.samples.RData")
load("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/pedigree/genotypedRelatives.RData")
row.names(grm) <- as.numeric(row.names(grm))-0.1
colnames(grm) <- row.names(grm)
names(grm) <- colnames(grm)

asib <- c()
bsib <- c()
for (i in seq_along(sibPairs)) {
    asib <- c(asib, as.character(sibPairs[[i]][1]))
    bsib <- c(bsib, as.character(sibPairs[[i]][2]))
}
aidx <- c()
bidx <- c()
for (i in 1:length(asib)) {
    aidx[[i]] <- which(row.names(grm) == asib[i])
    bidx[[i]] <- which(row.names(grm) == bsib[i])
}
sibValues <- c()
for (i in seq_along(aidx)) {
    sibValues[[i]] <- grm[aidx[i], bidx[i]]
}

ped.geno <- ped[ped$id %in% geno.samples$id,]
pedsib <- ped.geno[721:1830,]
pedsib$id <- pedsib$id - 0.1

sibgrm <- which(names(grm) %in% as.character(pedsib$id))
sibgrm <- grm[sibgrm,sibgrm]


plot(1:1229, sibgrm[1,])






library(lattice)
cols<- brewer.pal(11, 'PiYG')
test<-colorRampPalette(cols, alpha = TRUE)(50)
levelplot(sibgrm, col.regions=test, cuts=49, cex.axis=0.8,
          scales=list(x=list(rot=45)))




# sibships > 2
# 46043 46166 46203 48335 50277 50291 51967 51999 52129
# 63    90   102   137   265   268   436   442   464
# sibPairs[[510]] <- c(48355, 48417)
# sibPairs[[511]] <- c(48358, 48417)
# sibPairs[[512]] <- c(50150, 49217)
# sibPairs[[513]] <- c(50150, 49218)
# sibPairs[[514]] <- c(50150, 50150)
# sibPairs[[515]] <- c(50151, 49217)
# sibPairs[[516]] <- c(50151, 49218)
# sibPairs[[517]] <- c(49255, 49258)
# sibPairs[[518]] <- c(49255, 49234)
# sibPairs[[519]] <- c(49255, 49233)
# sibPairs[[520]] <- c(49233, 49258)
# sibPairs[[521]] <- c(49234, 49258)
# sibPairs[[522]] <- c(50147, 50189)
# sibPairs[[523]] <- c(50148, 50189)
# sibPairs[[524]] <- c(51304, 51295)
# sibPairs[[525]] <- c(51304, 51275)
# sibPairs[[526]] <- c(51313, 51216)
# sibPairs[[527]] <- c(51313, 51209)
# sibPairs[[528]] <- c(54315, 54304)
# sibPairs[[529]] <- c(54315, 54253)
# sibPairs[[530]] <- c(54315, 54244)
# sibPairs[[531]] <- c(54304, 54253)
# sibPairs[[532]] <- c(54304, 54244)
# sibPairs[[533]] <- c(54207, 54190)
# sibPairs[[534]] <- c(54207, 54187)
# sibPairs[[535]] <- c(54872, 54871)
# sibPairs[[536]] <- c(54872, 54817)
# sibPairs[[537]] <- c(54872, 54816)
# sibPairs[[538]] <- c(54871, 54817)
# sibPairs[[539]] <- c(54871, 54816)
#
# which(lapply(sibPairs, length) < 2)
# sibPairs <- sibPairs[-c(162,166,323,328,333)]


