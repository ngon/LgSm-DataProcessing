
## get covariate info for RNAseq samples

gex <- read.table("./brainSamples.ail.txt", sep="\t", as.is=T, header=T)[1:9]
phenos <- read.table("./dataFiles/phenotypes/allData.txt", sep="\t", header=T)

hip <- gex[gex$tissue == "hip",]
str <- gex[gex$tissue == "str",]
pfc <- gex[gex$tissue == "pfc",]

samples <- phenos[phenos$brain == "y",] # 313 mice whose brains were dissected
samples <- samples[,c("id", "gen", "dam", "sire", "dam.age", "sire.age", "wean.age", "mf.ratio", "sex", "batch", "rip.weight", "brain")]

sum(unique(str$id) %in% samples$id) # 231 unique out of 254 (23 dup)
sum(unique(hip$id) %in% samples$id) # 259 unique out of 280 (21 dup)
sum(unique(pfc$id) %in% samples$id) # 235 unique out of 255 (20 dup)

hipData <- merge(hip, samples, header=T, all.x=T)
strData <- merge(str, samples, header=T, all.x=T)
pfcData <- merge(pfc, samples, header=T, all.x=T)


################ IMPORTANT NOTE #####################
# The following 2 samples were sequenced but should be thrown out due to
# complications in phenotyping.

# 46025, the actual mouse with this id was discarded because she became
# pregnant during testing. the individual labeled as 46025 in the RNAseq data
# is likely to be a male with id 46130. during RNA extraction (pfc), it was noted # that the tube holding the brain sample was labeled 46130 on the cap, but 46025
# on the side. if 46025 is male (e.g. if he expresses SRY or other Y-chr mRNAs),
# keep him and relabel the mouse.id as 46130.

mouse46130 <- list(46130, "pfc", "pfc46130", pfcData[12,4:9], samples[48,2:13])
pfcData <- rbind(pfcData, unlist(mouse46130))
replace(pfcData[256,17], "M")
pfcData[256, 20] <- "y"

# 46348, also got pregnant and should not have sequenced at all. data for this
# mouse should be thrown out altogether.

gex.covariates <- rbind.data.frame(hipData, pfcData, strData)

# males are 0, females are 1
levels(gex.covariates$sex)[1:2] <- c(1,0)

save(gex.covariates, file="gex.covariates.Rdata")
write.table(gex.covariates, "./gex.covariates.txt", sep="\t", col.names=F,
            row.names=F, quote=F)


