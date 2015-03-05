# process phenotype for Natalia #

###### GO TO DIRECTORY
setwd("/group/palmer-lab/AIL/LgSm-DataProcessing")

###### CONSTANTS
MAF=0.05

#### Read genotype sample names
#### Change sample names -
#### 54109_11 changed to 54109 (double ear tag)
#### 51717-19 changed to 51717 (double ear tag)
#### 51348 changed to 51349 (inferred from pedigree checks) - possible typo - double check using sex(genotypes on X).
geno.samples <- read.table("/group/palmer-lab/AIL/GBS/dosage/genotyped.samples.txt", as.is=TRUE)
geno.samples[which(geno.samples$V1 == "54109_11"), 1] <- "54109"
geno.samples[which(geno.samples$V1 == "51717-19"), 1] <- "51717"
geno.samples[which(geno.samples$V1 == "51348"), 1]    <- "51349"
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
covar.names   <- read.table("cov.names.txt", skip=1, as.is=T)$V1
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

### covariates for each trait #### updated by Natalia 2/6/15
traitcovs <- vector("list", length=17)

names(traitcovs) <- c("cpp.diff", "act2.t", "act4.t", "wild.binary", "ppi6.logit", "ppi12.logit", 
                      "cpp8.1", "cpp8.t", "sens", "act1.t", "ppi3.logit",
                      "startle", "habituation", "act3.t", "act5.t", "act8.t","glucose")
                      
                      
traitcovs[["cpp.diff"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", 
                                "is.cpp.box7", "is.cpp.box11", "is.cpp.box12") 

traitcovs[["cpp8.1"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", 
                        "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12") 

traitcovs[["cpp8.t"]]  <- list("one", "sex", "is.gen53", "is.gen56", "is.batch2", "is.batch10", "is.batch13","is.batch15", "is.batch21", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12", "comerr8") 

traitcovs[["sens"]]    <- list("one", "sex", "is.gen52", "is.cpp.box7", "is.batch8") 

traitcovs[["act1.t"]]  <- list("one", "sex", "is.batch14", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "comerr1")

traitcovs[["act2.t"]]  <- list("one", "sex", "is.gen52","is.gen56", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box12", "is.batch8") 

traitcovs[["act4.t"]]  <- list("one", "sex", "is.cpp.box7", "is.cpp.box8", "is.cpp.box11") 

traitcovs[["wild.binary"]]    <- list("one", "sex")

traitcovs[["ppi3.logit"]] <- list("one", "sex", "is.ppi.box3", "ppi.weight", "is.batch4")
traitcovs[["ppi6.logit"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4", "ppi.weight")
traitcovs[["ppi12.logit"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4","ppi.weight", 
                            "is.batch3", "is.batch4", "is.batch7", "is.batch9")

traitcovs[["startle"]] <- list("one", "is.ppi.box3", "is.ppi.box4", "ppi.weight", "is.batch2", "is.batch16", "is.batch17") # shyam: is sex excluded for a reason - ## natalia: see comment under avg.ppi

traitcovs[["habituation"]] <- list("one", "ppi.weight") 

traitcovs[["act3.t"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55", "is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch21", "is.batch22", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12") 

traitcovs[["act5.t"]]  <- list("one", "sex", "is.gen51", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.cpp.box7", "is.cpp.box8", "is.batch7", "is.batch13", "is.batch15", "is.batch20", "is.batch21", "is.batch22") 

traitcovs[["act8.t"]]  <- list("one", "sex", "is.gen56", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "is.batch3", "is.batch7") 

traitcovs[["glucose"]] <- list("one", "sex", "glu.weight")


#################### Make script to run Gemma for chosen trait ##############
cmds <- c()
for (trait in names(traitcovs)) {
    index.pheno     <- which(pheno.names == trait)
    chosen.covars   <- covars[, unlist(traitcovs[[trait]])]
    write.table(chosen.covars, file=paste0("/group/palmer-lab/AIL/qtlmapping/covariates/", trait, ".covs"), sep="\t", quote=F, row.names=F, col.names=F)
    for (chrom in 1:19) {
        cmds <- c(cmds, paste0("gemma -g /group/palmer-lab/AIL/GBS/dosage/chr", chrom, ".filtered.dosage -p /group/palmer-lab/AIL/LgSm-DataProcessing/phenos.allgeno.txt -k /group/palmer-lab/AIL/qtlmapping/kinship/chrNot", chrom,".cXX.txt -a /group/palmer-lab/AIL/GBS/dosage/chr", chrom, ".filtered.snpinfo -c /group/palmer-lab/AIL/qtlmapping/covariates/", trait, ".covs -lmm 2 -maf ", MAF, " -o ", trait, ".chr", chrom, " -n ", index.pheno))
    }
}
write.table(cmds, file=paste0("/group/palmer-lab/AIL/code/gemma.alltraits.cmds"), row.names=F, col.names=F, quote=F)

