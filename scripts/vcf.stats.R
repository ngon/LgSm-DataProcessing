#### 10 Jul 2015 --------------------------------------------------------------
# Summary stats on variants called from known SNPs in the Cheverud data set.
# Stat files were generated using vcftools commands on the files in /GBS/vars.

### MAKE COMMAND FILES FOR VCFTOOLS
# chrom <- paste0("chr", 1:19)
# ld.cmds <- c()
# for (c in chr){
#     ld.cmds <- c(ld.cmds, paste0("vcftools --vcf ail.", c, ".known.vcf --out ./vcf.stats/", c,
#                                  " --geno-r2 --ld-window-bp 10000000 --min-r2 0"))
# }
# cmds <- c()
# for (c in chrom){
#     cmds <- c(cmds, paste0("vcftools --vcf ail.", c, ".known.vcf --out ./vcf.stats/", c,
#                            "--freq --counts --depth --site-depth --site-mean-depth ",
#                            "--geno-depth --site-quality --het --hardy --missing --SNPdensity ",
#                            "--TsTv-by-count --TsTv-by-qual --filtered-sites --singletons"))
# }

#### PREPARE DATA -------------------------------------------------------------

### SUMMARY STATISTICS PER INDIVIDUAL ###

# individual read depth
# replace IDs with ids in geno.samples
idepth <- read.table("./dataFiles/chr19.idepth", header=T, sep="\t", as.is=T)
ids <- read.table("./genotyped.samples.txt", sep="\t", header=F)[1]
names(ids) <- "INDV"
idepth$INDV <- ids$INDV

# missingness per individual
imiss <- read.table("./dataFiles/chr19.imiss", sep="\t", header=T, as.is=T)
imiss$INDV <- ids$INDV

# heterozygosity per individual
# 15 individuals were excluded, so we can't replace the first col with geno.
# samples. we also have to replace ids with an underscore and ids that were
# mislabeled (typos).
het <- read.table("./dataFiles/chr19.het", header=T, sep="\t", as.is=T)
id <- unlist(strsplit(het$INDV, "[.]"))
id <- id[which(seq_along(id) %% 2 != 0)]
id[1383] <- "51717" # underscore
id[1591] <- "54109" # underscore
id[1815] <- "51801" # typo
id[1301] <- "51349" # typo
id <- as.integer(id)
het$INDV <- id
names(het)[4] <- "N_SITES_HET"

# get IDs of mice that were excluded
excludedMice <- idepth$INDV[which(!idepth$INDV %in% het$INDV)]

# merge into one data frame and clean up
perMouse <- merge(idepth, imiss, all=T)
perMouse <- merge(perMouse, het, by="INDV", all.x=T)
rm(het, idepth, imiss)

### SUMMARY STATISTICS PER LOCUS ###

# get true singletons ("S") and private doubletons ("D"). doubletons are SNPs
# where the minor allele only occurs in 1 mouse, and the mouse is homozygous.
singletons <- read.table("./dataFiles/chr19.singletons", header=T, sep="\t", as.is=T)
id <- unlist(strsplit(singletons$INDV, "[.]"))
id <- id[which(seq_along(id) %% 2 != 0)]
which(id == "57801")
which(id == "51348")
which(id == "51717_19")
which(id == "54109_11")
id[959] <- "54109"
id <- as.integer(id)
singletons$INDV <- id
names(singletons)[3] <- "SING.DOUB"

# read depth
ldepth <- read.table("./dataFiles/chr19.ldepth", header=T, sep="\t", as.is=T)
ldepthMean <- read.table("./dataFiles/chr19.ldepth.mean", header=T, sep="\t", as.is=T)
identical(ldepth$POS, ldepthMean$POS)
depth <- cbind(ldepth, ldepthMean$MEAN_DEPTH, ldepthMean$VAR_DEPTH)
names(depth)[5:6] <- c("MEAN_DEPTH", "VAR_DEPTH")

# hardy weinberg equilibrium
hwe <- read.table("./dataFiles/chr19.hwe", header=T, sep="\t", as.is=T)
obs <- unlist(strsplit(hwe$OBS.HOM1.HET.HOM2., "[/]"))
O.HOM1 <- obs[seq(1, length(obs), by=3)]
O.HET <- obs[(2, length(obs), by=3)]
O.HOM2 <- obs[seq(3, length(obs), by=3)]
exp <- unlist(strsplit(hwe$E.HOM1.HET.HOM2., "[/]"))
E.HOM1 <- exp[seq(1, length(exp), by=3)]
E.HET <- exp[seq(2, length(exp), by=3)]
E.HOM2 <- exp[seq(3, length(exp), by=3)]
hardy <- cbind(hwe[2], O.HOM1, O.HET, O.HOM2, E.HOM1, E.HET, E.HOM2, hwe[5:6])

# merge and clean
perLocus <- merge(depth, hardy, by="POS", all.x=T)
rm(exp, obs, O.HOM1, O.HOM2, O.HET, E.HET, E.HOM1, E.HOM2, hardy, hwe, depth,
   ldepthMean, ldepth, id)

# missingness and quality per site (QUAL column of VCF file)
lqual <- read.table("./dataFiles/chr19.lqual", header=T, sep="\t", as.is=T)
lmiss <- read.table("./dataFiles/chr19.lmiss", header=T, sep="\t", as.is=T)
identical(lqual$POS, lmiss$POS)
lmiss <- cbind(lmiss, lqual$QUAL)
names(lmiss)[7] <- "QUAL"

# merge and clean
perLocus <- merge(perLocus, lmiss, all.x=T, all.y=T)
rm(lmiss, lqual)

# allele frequency
frq <- read.table("./dataFiles/chr19.frq", header=T, sep="\t", as.is=T)
# replace blank values in ALLELE.FREQ3 with NA
all3 <- c()
for (f in seq_along(frq$ALLELE.FREQ3)) {
    if (frq$ALLELE.FREQ3[f] == "") { all3 <- c(all3, "NA:NA") }
    if (frq$ALLELE.FREQ3[f] != "") { all3 <- c(all3, frq$ALLELE.FREQ3[f]) }
}

# there are 2 rows in frq for which ALLELE.FREQ1 and ALLELE.FREQ2 values are
# blank. find these and change them to NA.
# test<- grep(pattern="[:]", x=frq$ALLELE.FREQ1)
# a<-frq[!seq_along(frq$ALLELE.FREQ1) %in% test,]
# test2<- grep(pattern="[:]", x=frq$ALLELE.FREQ2)
# b<-frq[!seq_along(frq$ALLELE.FREQ2) %in% test,]
frq$ALLELE.FREQ1[438] <- "NA:NA"
frq$ALLELE.FREQ1[3318] <- "NA:NA"
frq$ALLELE.FREQ2[438] <- "NA:NA"
frq$ALLELE.FREQ2[3318] <- "NA:NA"

a1 <- unlist(strsplit(frq$ALLELE.FREQ1, "[:]"))
a2 <- unlist(strsplit(frq$ALLELE.FREQ2, "[:]"))
a3 <- unlist(strsplit(all3, "[:]"))
ALLELE1 <- a1[seq(1, length(a1), by=2)]
ALLELE2 <- a2[seq(1, length(a2), by=2)]
ALLELE3 <- a3[seq(1, length(a3), by=2)]
FREQ1 <- a1[seq(2, length(a1), by=2)]
FREQ2 <- a2[seq(2, length(a2), by=2)]
FREQ3 <- a3[seq(2, length(a3), by=2)]

frq <- cbind(frq[2:4], ALLELE1, ALLELE2, ALLELE3, FREQ1, FREQ2, FREQ3)

# allele count
frqCount <- read.table("./dataFiles/chr19.frq.count", header=T, sep="\t", as.is=T)
all3 <- c()
for (f in seq_along(frqCount$ALLELE.COUNT3)) {
    if (frqCount$ALLELE.COUNT3[f] == "") { all3 <- c(all3, "NA:NA") }
    if (frqCount$ALLELE.COUNT3[f] != "") { all3 <- c(all3, frqCount$ALLELE.COUNT3[f]) }
}
a1 <- unlist(strsplit(frqCount$ALLELE.COUNT1, "[:]"))
a2 <- unlist(strsplit(frqCount$ALLELE.COUNT2, "[:]"))
a3 <- unlist(strsplit(all3, "[:]"))
COUNT1 <- a1[seq(2, length(a1), by=2)]
COUNT2 <- a2[seq(2, length(a2), by=2)]
COUNT3 <- a3[seq(2, length(a3), by=2)]

frqCount <- cbind(frqCount[2], COUNT1, COUNT2, COUNT3)

# merge and cleanup
freq <- merge(frq, frqCount, by="POS", all.x=T)
perLocus <- merge(perLocus,freq, by="POS", all.x=T, all.y=T)
perLocus$CHR <- NULL
rm(a1, a2, a3, all3, ALLELE1, ALLELE2, ALLELE3, COUNT1, COUNT2, COUNT3, FREQ1,
   FREQ2, FREQ3, f, freq, frq, frqCount)

# save the organized data frames.
save(excludedMice, perMouse, singletons, perLocus, file="./dataFiles/vcfStats.RData")


####### EXPLORATORY ANALYSIS ---------------------------------------------------

### hwe is obv way crazy without filtering on missingness.
pvals.o <- -log10(sort(perLocus$P, decreasing=F))
pvals.o <- pvals.o[2021:length(pvals.o)]
pvals.e <- -log10(1:length(pvals.o)/length(pvals.o))
library(ggplot2)
plot(pvals.e,pvals.o,pch=19,cex=0.25,
         main="Expected vs Observed HWE p-values",
         xlab=expression(Expected~~-log[10](italic(p))),
         ylab=expression(Observed~~-log[10](italic(p))),
         xlim=c(0,max(pvals.o)),
         ylim=c(0,max(pvals.o)))
lines(pvals.e,pvals.e,col="red")


23^2

##### NOT USING THESE FILES AT THE MOMENT ##### -------------------------------
# load kept and removed sites: useless here bc i didn't specify any filters.
# kept.sites <- read.table("./dataFiles/chr19.kept.sites", header=T, sep="\t", as.is=T)
# removed.sites <- read.table("./dataFiles/chr19.removed.sites", header=T, sep="\t", as.is=T)
# there are 7991 SNPs on chr 19

# transition to transversion ratio, by alt.allele count and quality
#tstv.qual <- read.table("./dataFiles/chr19.TsTv.qual", header=T, sep="\t", as.is=T)
#tstv.count <- read.table("./dataFiles/chr19.TsTv.count", header=T, sep="\t",as.is=T)

# genotype depth with missing entries coded as -1. snp x indiv matrix.
# gdepth <- read.table("./dataFiles/chr19.gdepth", header=T, sep="\t", as.is=T,
#                      na.strings = "-1")
# names(gdepth)[3:1832] <- ids$INDV
# row.names(gdepth) <- gdepth$POS
# gdepth <- gdepth[-c(1,2)]


