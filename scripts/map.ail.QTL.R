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

### covariates for each trait #### updated by Natalia 2/6/15
traitcovs <- vector("list", length=17)

names(traitcovs) <- c("cpp.diff", "cpp8.1", "cpp8.t", "sens",
                      "act1.t", "act2.t", "act4.t", "wild", "ppi3", "ppi6", "ppi12",
                      "startle", "habituation", "act3.t", "act5.t", "act8.t","glucose")
                      
                      
traitcovs[["cpp.diff"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", 
                                "is.cpp.box7", "is.cpp.box11", "is.cpp.box12") 

traitcovs[["cpp8.1"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", 
                        "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12") 

traitcovs[["cpp8.t"]]  <- list("one", "sex", "is.gen53", "is.gen56", "is.batch2", "is.batch10", "is.batch13","is.batch15", "is.batch21", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12") 

traitcovs[["sens"]]    <- list("one", "sex", "is.gen52", "is.cpp.box7", "is.batch8") 

traitcovs[["act1.t"]]  <- list("one", "sex", "is.batch14", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12")

traitcovs[["act2.t"]]  <- list("one", "sex", "is.gen52","is.gen56", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box12", "is.batch8") 

traitcovs[["act4.t"]]  <- list("one", "sex", "is.cpp.box7", "is.cpp.box8", "is.cpp.box11") 

traitcovs[["wild"]]    <- list("one", "sex")

traitcovs[["ppi3"]] <- list("one", "sex", "is.ppi.box3", "ppi.weight", "is.batch4")
traitcovs[["ppi6"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4", "ppi.weight")
traitcovs[["ppi12"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4","ppi.weight", 
                            "is.batch3", "is.batch4", "is.batch7", "is.batch9")

traitcovs[["startle"]] <- list("one", "is.ppi.box3", "is.ppi.box4", "ppi.weight", "is.batch2", "is.batch16", "is.batch17") # shyam: is sex excluded for a reason - ## natalia: see comment under avg.ppi

traitcovs[["habituation"]] <- list("one", "ppi.weight") 

traitcovs[["act3.t"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55", "is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch21", "is.batch22", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12") 

traitcovs[["act5.t"]]  <- list("one", "sex", "is.gen51", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.cpp.box7", "is.cpp.box8", "is.batch7", "is.batch13", "is.batch15", "is.batch20", "is.batch21", "is.batch22") 

traitcovs[["act8.t"]]  <- list("one", "sex", "is.gen56", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "is.batch3", "is.batch7") 

traitcovs[["glucose"]] <- list("one", "sex", "glu.weight")


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

