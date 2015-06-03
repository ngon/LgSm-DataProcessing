##### MAP QTL IN F50-56 AIL #####
# This file uses permuted phenotype data (5 permutations for each of 14 traits).

##### GO TO DIRECTORY
setwd("/group/palmer-lab/AIL/LgSm-DataProcessing")

##### CONSTANTS
MAF=0.05

#### READ GENOTYPED SAMPLES AND CHANGE SOME NAMES  ----------------------
# 54109_11 changed to 54109 (double ear tag)
# 51717-19 changed to 51717 (double ear tag)
# 51348 changed to 51349 (inferred from ped checks - typo? - double check using sex(genotypes on X).
geno.samples <- read.table("/group/palmer-lab/AIL/GBS/dosage/genotyped.samples.txt", as.is=TRUE)
geno.samples[which(geno.samples$V1 == "54109_11"), 1] <- "54109"
geno.samples[which(geno.samples$V1 == "51717-19"), 1] <- "51717"
geno.samples[which(geno.samples$V1 == "51348"), 1]    <- "51349"
names(geno.samples) <- c("id")

# READ TAB DELIMITED FILE WITH PERMUTED PHENOTYPES
# Sex is an indicator variable; M=0 and F=1
pheno <- read.table("./permutedPhenotypes.txt", sep="\t", header=T, as.is=T)
#covars <- read.table("./covariates.orm.txt", sep="\t", header=T, as.is=T)
pheno.names <- names(pheno)
ids<- read.table("./dataFiles/", header=T, as.is=T)[2]
names(ids)[1] <- "id"
pheno <- cbind(ids[1], pheno)

#### EXTRACT DATA FOR GENOTYPED SAMPLES ------------------------------------
### Generate phenotype data for all samples with genotypes
permPheno.allgeno <- merge(geno.samples, pheno, all.x=TRUE)
write.table(permPheno.allgeno, file="permPhenos.allgeno.txt",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


### Make covariate file with all ids in genotype file
traits.to.perm <- c("act1.t", "act2.t", "act4.t", "act5.t", "act8.t", "cpp.diff",
                    "sc8.t", "sc1.t", "startle", "wild.binary", "ppi6.logit",
                    "ppi12.logit", "tail","glucose")

get.filenames <- function(trait){
filenames <- c()
for (i in 1:5){
filenames[i] <- paste0("/group/palmer-lab/AIL/qtlmapping/covariates/permCovs/perm",
                    i, ".", trait, ".cov.txt")
}
return(filenames)
}
files <- lapply(traits.to.perm, get.filenames)
files <- unlist(files)

merge.covars <- function(covFile, genoFile, ids){
    covariate <- read.table(covFile, as.is=T, header=T)
    covariate <- cbind(ids, covariate)
    covariate <- merge(genoFile, covariate, header=T, all.x=TRUE)
    covariate$one <- 1
    write.table(covariate[-1], file=paste0(covFile,"2"), sep="\t", row.names=F, col.names=F, quote=F)
}

lapply(files, merge.covars, genoFile=geno.samples, ids=ids)


### MAKE SCRIPT TO RUN GEMMA FOR CHOSEN TRAIT ----------------------------------------------

cmds <- c()
for (pheno in pheno.names) {

    index.pheno     <- which(pheno.names == pheno)
    chosen.covars   <- paste0("/group/palmer-lab/AIL/qtlmapping/covariates/permCovs/",
                              pheno, ".cov.txt2")


    for (chrom in 1:19) {
        cmds <- c(cmds, paste0("gemma -g /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr", chrom,
                               ".filtered.dosage -p /group/palmer-lab/AIL/LgSm-DataProcessing/permPhenos.allgeno.txt -k /group/palmer-lab/AIL/qtlmapping/kinship/onlyEmpirical/chrNot",
                               chrom,".cXX.txt -a /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr",
                               chrom, ".filtered.snpinfo -c ", chosen.covars,
                               " -lmm 2 -maf ", MAF, " -o ", pheno, ".chr", chrom, " -n ", index.pheno))
}
write.table(cmds, file=paste0("/group/palmer-lab/AIL/code/gemma.permPheno.emp.cmds"),
            row.names=F, col.names=F, quote=F)
}
