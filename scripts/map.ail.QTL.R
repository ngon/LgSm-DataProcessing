# process phenotype for Natalia #

###### CONSTANTS
MAF=0.05

#### Read genotype sample names
geno.samples <- read.table("genotyped.samples.txt", as.is=TRUE)
names(geno.samples) <- c("id")

#### Read phenotype file
pheno     <- read.table("pheno.noHeader.txt", as.is=TRUE)
pheno.names <- read.table("pheno.names.txt", skip=1, as.is=T)[,2]
names(pheno) <- pheno.names

#### Generate the phenotype data for all samples with genotypes
pheno.allgeno <- merge(geno.samples, pheno, all.x=TRUE)
write.table(pheno.allgeno, file="phenos.allgeno.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

#### Read master covariate file
covars        <- read.table("cov.noHeader.txt", sep="\t", as.is=T)
covar.names   <- read.table("cov.names.txt", skip=1, as.is=T)[,2]
names(covars) <- covar.names

#### Change the sex covariate to an indicator variable
covars$sex    <- replace(covars$sex, covars$sex == "M", "0")
covars$sex    <- replace(covars$sex, covars$sex == "F", "1")
covars$sex    <- as.numeric(covars$sex)

#### Change the batch and ppi.box covariates to set of indicator variables for gemma
#### Convert generation to factor as well.
batches       <- unique(covars$batch)
#batches       <- batches[-1]
for (batch in batches) {
  covars[,paste0("is.batch", batch)]  <- as.integer(covars$batch == batch)
}
boxes         <- unique(covars$ppi.box[which(!is.na(covars$ppi.box))])
#boxes         <- boxes[-1]
for (box in boxes) {
  covars[,paste0("is.ppi.box", box)] <- as.integer(covars$ppi.box == box)
}
boxes         <- unique(covars$cpp.box[which(!is.na(covars$cpp.box))])
#boxes         <- boxes[-1]
for (box in boxes) {
  covars[,paste0("is.cpp.box", box)] <- as.integer(covars$cpp.box == box)
}
gens          <- unique(covars$gen)
#gens          <- gens[-1]
for (gen in gens) {
    covars[,paste0("is.gen", gen)] <- as.integer(covars$gen == gen)
}

#### Make covariate file with all ids in gentype file
covars <- merge(geno.samples, covars, all.x=TRUE)
covars$one <- 1

####### Code for anova testing batch having a significant explantory power on outcome
####### and batch 13 as being different from rest of batches in terms of outcome
#cpp8.t <- pheno$cpp8.t
#batch <- as.factor(covars$batch)
#anova.cpp <- anova(lm(cpp8.t~batch))
#anova.cpp <- anova(lm(cpp8.t~covars$is.batch13))
#anova.cpp <- anova(lm(cpp8.t~as.factor(covars$cpp.box)))

### covariates for each trait
traitcovs <- vector("list", length=12)
names(traitcovs) <- c("cpp8.t", "sens", 
                      "act1.t", "act2.t", "act3.t", "act4.t", "act5.t", "act8.t",
                      "avg.ppi", "startle", "glucose", "wild")

traitcovs[["cpp8.t"]]  <- list("one", "sex", "is.gen53", "is.batch2", "is.batch10", "is.cpp.box1", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12") # gen53,56, batch 2,10,13,15,21, all except box 5 - "is.gen56",  "is.batch13", "is.batch15", "is.batch21", - removed due to no genotyped samples from these batches/gens
traitcovs[["sens"]]    <- list("one", "sex", "is.gen52", "is.cpp.box7", "is.batch8") # gen52, box7, batch 8
traitcovs[["act1.t"]]  <- list("one", "sex", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12")#, "is.batch14") # box3-9,11,12, batch 14 - removed batch 14 for now due to no samples from batch 14 in genotyped samples
traitcovs[["act2.t"]]  <- list("one", "sex", "is.gen52", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box12", "is.batch8") # gen52,56, box5,7,8,10,12, batch 8 - removed is.gen56 as no genotyped samples
traitcovs[["act3.t"]]  <- list("one", "sex", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12") # box7,8,10,11,12 - shyam removed gen and batch since which gens and batches are correlated with outcome was not mentioned in comments by natalia
traitcovs[["act4.t"]]  <- list("one", "sex", "is.cpp.box7", "is.cpp.box8", "is.cpp.box11") # box7,8,11
traitcovs[["act5.t"]]  <- list("one", "sex", "is.gen51", "is.gen53", "is.cpp.box7", "is.cpp.box8", "is.batch7") # box 7,8, gen 51,53-56, batch 7,13,15,20-22 -  "is.gen54", "is.gen55", "is.gen56", "is.batch13", "is.batch15", "is.batch20", "is.batch21", "is.batch22" - removed due to no samples from this batch / gen in genotyped data
traitcovs[["act8.t"]]  <- list("one", "sex", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "is.batch3", "is.batch7") # box 3,4,6-9,11,12, gen56, batch 3,7 -  "is.gen56", removed due to no genotyped samples from gen 56
traitcovs[["avg.ppi"]] <- list("one", "is.ppi.box3", "is.ppi.box4", "ppi.weight", "is.batch4") #box 3, box4, batch4 - is sex excluded for a reason
traitcovs[["startle"]] <- list("one", "is.ppi.box3", "is.ppi.box4", "ppi.weight", "is.batch2") # box 3,4, batch 2,16,17 - is sex excluded for a reason - , "is.batch16", "is.batch17" - removed as no genotyped samples from these batches
traitcovs[["glucose"]] <- list("one", "sex", "glu.weight")
traitcovs[["wild"]]    <- list("one", "sex")

#################### Make script to run Gemma for chosen trait ##############
for (trait in names(traitcovs)) {
    index.pheno     <- which(pheno.names == trait)
    chosen.covars   <- covars[, unlist(traitcovs[[trait]])]
    write.table(chosen.covars, file=paste0("covariates/", trait, ".covs"), sep="\t", quote=F, row.names=F, col.names=F)
    cmds <- ("#! /bin/bash\n#$ -cwd\n#$ -j y\n")
    for (chrom in 1:19) {
        cmds <- c(cmds, paste0("/home/shyamg/bin/gemma -g genotypes/ail.chr", chrom, ".filtered.dosage -p phenos.allgeno.txt -k kinship/notChr", chrom,".cXX.txt -a snpinfo/ail.chr", chrom, ".snpinfo -c covariates/", trait, ".covs -lmm 2 -maf ", MAF, " -o ", trait, ".chr", chrom, " -n ", index.pheno))
    }
    write.table(cmds, file=paste0("scripts/gemma.", trait, ".sh"), row.names=F, col.names=F, quote=F)
}

#group.name <- "cpp8"
#chosen.pheno <- c("cpp8.t")
#chosen.covars <- c("one", "sex", paste0("is.cpp.box", c(1:12)[-5]))
#index.pheno <- sapply(chosen.pheno, function(x) { which(pheno.names == x)} )
## Make covariate file.
#chosen.covars <- sapply(chosen.covars, function(x) { which(names(covars) == x)} )
#chosen.covars <- covars[, chosen.covars]
#write.table(chosen.covars, file=paste0("covariates/", group.name, ".covs"), sep="\t", quote=F, row.names=F, col.names=F)

## Run gemma chromosome wise
#cmds <- c("#! /bin/bash\n#$ -cwd\n#$ -j y\n")
#for (index in 1:length(index.pheno)) {
#    for (chrom in 1:19) {
#        cmds <- c(cmds, paste0("/home/shyamg/bin/gemma -g genotypes/ail.chr", chrom, ".dosage -p phenosAIL.allgeno.txt -k kinship/notChr", chrom,".cXX.txt -a snpinfo/ail.chr", chrom, ".snpinfo -c covariates/", group.name, ".covs -lmm 2 #-maf ", MAF, " -o ", chosen.pheno[index], ".chr", chrom, ".out -n ", index.pheno[index]))
#    }
#    write.table(cmds, file=paste0("scripts/", group.name, ".", chosen.pheno[index], ".sh"), row.names=F, col.names=F, quote=F)
#}

