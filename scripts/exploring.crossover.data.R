
## EXPLORING CROSSOVER DATA ##

# get ped info for genotyped samples
ped <- read.table("./pedigree/pedforQTLRel.txt", header=T)
load("./pedigree/geno.samples.RData")

names(ped)[5] <- "gen"
ped <- ped[ped$id %in% geno.samples$id,]
rm(geno.samples)

# get crossover data
source("./scripts/hap.to.color.R")
load("./dataFiles/chr19.hap.RData")

hapdata <- hap.to.color(hap=chr19.hap)
names(hapdata) <- c("breaks", "col.breaks")
rm(chr19.hap)

lbreaks<- hapdata[["breaks"]][["left"]]
rbreaks <- hapdata[["breaks"]][["right"]]

# crossover stats
llens <- unlist(lapply(lbreaks, length))
rlens <- unlist(lapply(rbreaks, length))
mean(llens) # 3464.683
mean(rlens) #3361.686
quantile(llens, probs=seq(0, 1, 0.1))
quantile(rlens, probs=seq(0, 1, 0.1))

# get crossover stats for each generation
# male and female crossover stats by generation
# compare crossovers within families


