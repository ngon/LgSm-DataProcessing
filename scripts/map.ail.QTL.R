##### MAP QTL IN F50-56 AIL #####
## This file uses phenotype data with outliers removed.
## For details on how this was done, see ./dataFiles/inspectOutliers/README.txt and
## look at the script ./cookieTime/inspectOutliers.R


##### GO TO DIRECTORY
setwd("/group/palmer-lab/AIL/LgSm-DataProcessing")

##### CONSTANTS
MAF=0.05


#### READ GENO, PHENO & COV FILES -----------------------------------------------

#### READ GENOTYPED SAMPLES AND CHANGE SOME NAMES  ----------------------
# 54109_11 changed to 54109 (double ear tag)
# 51717-19 changed to 51717 (double ear tag)
# 51348 changed to 51349 (inferred from ped checks - typo? - double check using sex(genotypes on X).
geno.samples <- read.table("/group/palmer-lab/AIL/GBS/dosage/genotyped.samples.txt", as.is=TRUE)
geno.samples[which(geno.samples$V1 == "54109_11"), 1] <- "54109"
geno.samples[which(geno.samples$V1 == "51717-19"), 1] <- "51717"
geno.samples[which(geno.samples$V1 == "51348"), 1]    <- "51349"
names(geno.samples) <- c("id")


# READ TAB DELIMITED PHENO AND COV FILES ** WITH OUTLIERS REMOVED **
# Sex is an indicator variable; M=0 and F=1
pheno <- read.table("./phenotypes.orm.txt", sep="\t", header=T, as.is=T)
covars <- read.table("./covariates.orm.txt", sep="\t", header=T, as.is=T)
pheno.names <- names(pheno)

#### EXTRACT DATA FOR GENOTYPED SAMPLES ------------------------------------

### Generate phenotype data for all samples with genotypes
pheno.allgeno <- merge(geno.samples, pheno, all.x=TRUE)
write.table(pheno.allgeno, file="phenos.allgeno.txt",
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

### Make covariate file with all ids in genotype file
covars <- merge(geno.samples, covars, all.x=TRUE)
#covars$one <- 1

#### DEFINE COVARIATES FOR EACH TRAIT ------------------------------------

# GET TRAIT NAMES
# ppi traits are logit-transformed.
# wild, is.coatA, is.coatB, and is.coatW are binary.

traitcovs <- vector("list", length=84)

names(traitcovs) <- c("ppi3.logit", "ppi6.logit", "ppi12.logit", "wild.binary", "tail","glucose",
                      "sc1.t", "sc1.1", "sc1.2", "sc1.3", "sc1.4", "sc1.5", "sc1.6",
                      "sc8.t", "sc8.1", "sc8.2", "sc8.3", "sc8.4", "sc8.5", "sc8.6",
                      "cpp1.t", "cpp1.1", "cpp1.2", "cpp1.3", "cpp1.4", "cpp1.5", "cpp1.6",
                      "cpp8.t", "cpp8.1", "cpp8.2", "cpp8.3", "cpp8.4", "cpp8.5", "cpp8.6",
                      "act1.t", "act1.1", "act1.2", "act1.3", "act1.4", "act1.5", "act1.6",
                      "act2.t", "act2.1", "act2.2", "act2.3", "act2.4", "act2.5", "act2.6",
                      "act3.t", "act3.1", "act3.2", "act3.3", "act3.4", "act3.5", "act3.6",
                      "act4.t", "act4.1", "act4.2", "act4.3", "act4.4", "act4.5", "act4.6",
                      "act5.t", "act5.1", "act5.2", "act5.3", "act5.4", "act5.5", "act5.6",
                      "act8.t", "act8.1", "act8.2", "act8.3", "act8.4", "act8.5", "act8.6",
                      "habituation", "startle", "sens", "cpp.diff","cpp.diff.p",
                      "is.coatA", "is.coatB", "is.coatW"
)

# PPI, STARTLE AND HABITUATION
traitcovs[["ppi3.logit"]] <- list("one", "sex", "is.ppi.box3", "ppi.weight", "is.batch4")
traitcovs[["ppi6.logit"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4", "ppi.weight")
traitcovs[["ppi12.logit"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4","ppi.weight",
                                   "is.batch3", "is.batch4", "is.batch7", "is.batch9")
traitcovs[["habituation"]] <- list("one", "sex", "ppi.weight")
traitcovs[["startle"]] <- list("one", "sex", "is.ppi.box3", "is.ppi.box4", "ppi.weight",
                               "is.batch2", "is.batch16", "is.batch17")


# SIDE CHANGES DAY 1 (total and six 5-min bins)
traitcovs[["sc1.t"]] <-list("one", "sex", "is.gen51", "is.gen53", "is.gen56", "is.gen54", "is.gen55","is.batch2", "is.batch4", "is.batch6", "is.batch12", "is.batch13", "is.batch21", "is.cpp.box2", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box12")
traitcovs[["sc1.1"]] <-list("one", "sex", "is.gen51", "is.gen56", "is.batch2", "is.batch4", "is.batch6", "is.batch12", "is.batch13", "is.batch20", "is.batch21","is.cpp.box6", "is.cpp.box11")
traitcovs[["sc1.2"]] <-list("one", "sex", "is.gen51", "is.gen53", "is.gen56","is.batch2", "is.batch4", "is.batch6", "is.batch12", "is.batch13", "is.batch21", "is.batch22", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box12")
traitcovs[["sc1.3"]] <-list("one", "sex", "is.gen51", "is.gen53", "is.gen56", "is.batch2", "is.batch4", "is.batch6", "is.batch13", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box12")
traitcovs[["sc1.4"]] <-list("one", "sex", "is.gen51", "is.gen53", "is.gen56", "is.gen54","is.batch4", "is.batch6", "is.batch12", "is.batch13", "is.batch21", "is.cpp.box2", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7","is.cpp.box8", "is.cpp.box10", "is.cpp.box12")
traitcovs[["sc1.5"]] <-list("one", "sex", "is.gen51", "is.gen53", "is.batch4", "is.batch6", "is.batch12", "is.batch13", "is.batch21", "is.cpp.box6", "is.cpp.box11")
traitcovs[["sc1.6"]] <-list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen56", "is.gen54", "is.gen55","is.batch4", "is.batch6","is.batch9", "is.batch13", "is.batch19", "is.batch21", "is.cpp.box6", "is.cpp.box11", "is.cpp.box12")


# SIDE CHANGES DAY 8 (total and six 5-min bins)
traitcovs[["sc8.t"]] <-list("one", "sex", "is.gen54", "is.gen55", "is.gen56", "is.batch2","is.batch3", "is.batch6", "is.batch7", "is.batch10","is.batch15", "is.batch20", "is.cpp.box6")
traitcovs[["sc8.1"]] <-list("one", "sex", "is.gen51", "is.batch3", "is.batch6", "is.batch19", "is.batch20", "is.cpp.box8")
traitcovs[["sc8.2"]] <-list("one", "sex", "is.gen51", "is.gen53","is.gen55", "is.gen56", "is.batch2","is.batch3", "is.batch7", "is.batch11","is.batch14","is.batch15", "is.cpp.box2", "is.cpp.box4", "is.cpp.box6", "is.cpp.box11")
traitcovs[["sc8.3"]] <-list("one", "sex", "is.gen51","is.batch20", "is.cpp.box2", "is.cpp.box6", "is.cpp.box12")
traitcovs[["sc8.4"]] <-list("one", "sex", "is.gen51", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.batch6", "is.batch15", "is.cpp.box11")
traitcovs[["sc8.5"]] <-list("one", "sex", "is.gen51", "is.gen54", "is.gen55", "is.gen56", "is.batch3", "is.batch15", "is.batch20", "is.cpp.box6")
traitcovs[["sc8.6"]] <-list("one", "sex", "is.gen51", "is.gen53", "is.gen55", "is.batch2","is.batch3", "is.batch6", "is.batch8", "is.cpp.box4", "is.cpp.box6", "is.cpp.box12")


# CPP DAY 1 - INITIAL PREFERENCE (total and six 5-min bins)
traitcovs[["cpp1.t"]]  <- list("one", "sex", "is.gen53", "is.gen54", "is.gen55","is.gen56", "is.batch11", "is.batch13", "is.batch14","is.batch17", "is.batch19","is.batch21", "is.batch22", "is.cpp.box3", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11")
traitcovs[["cpp1.1"]]  <- list("one", "sex", "is.gen55","is.gen56", "is.batch2", "is.batch7", "is.batch14","is.batch17", "is.batch19","is.batch21", "is.batch22", "is.cpp.box5", "is.cpp.box8")
traitcovs[["cpp1.2"]]  <- list("one", "sex", "is.gen53", "is.gen54", "is.gen56", "is.batch2", "is.batch6", "is.batch11","is.batch13", "is.batch14","is.batch20", "is.batch22", "is.cpp.box5")
traitcovs[["cpp1.3"]]  <- list("one", "sex", "is.batch7", "is.batch13", "is.batch18", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11")
traitcovs[["cpp1.4"]]  <- list("one", "sex", "is.gen56", "is.batch13", "is.batch22", "is.cpp.box4", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10")
traitcovs[["cpp1.5"]]  <- list("one", "sex", "is.gen51", "is.gen54", "is.gen55","is.gen56", "is.batch2", "is.batch4", "is.batch13","is.batch14", "is.batch17","is.batch18", "is.batch20", "is.batch22", "is.cpp.box11")
traitcovs[["cpp1.6"]]  <- list("one", "sex", "is.gen53", "is.gen54", "is.gen55","is.gen56", "is.batch11", "is.batch13","is.batch17", "is.batch20","is.batch21", "is.batch22", "is.cpp.box3", "is.cpp.box8", "is.cpp.box10")


# CPP DAY 8 - CONDITIONED PREFERENCE (total and six 5-min bins)
traitcovs[["cpp8.t"]]  <- list("one", "sex", "is.gen53", "is.gen56", "is.batch2", "is.batch10", "is.batch13","is.batch15", "is.batch21", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12", "comerr8")
traitcovs[["cpp8.1"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["cpp8.2"]]  <- list("one", "sex", "is.gen53", "is.batch3", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box8", "is.cpp.box11", "is.cpp.box12")
traitcovs[["cpp8.3"]]  <- list("one", "sex","is.batch15", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box8", "is.cpp.box11", "is.cpp.box12", "comerr8")
traitcovs[["cpp8.4"]]  <- list("one", "sex","is.gen53", "is.batch2","is.batch10","is.batch15","is.batch21", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12")
traitcovs[["cpp8.5"]]  <- list("one", "sex", "is.batch2", "is.batch5", "is.batch10","is.batch13", "is.batch14", "is.batch17", "is.batch18", "is.batch19", "is.batch21",  "is.gen51", "is.gen56", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box7","is.cpp.box8", "is.cpp.box9","is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["cpp8.6"]]  <- list("one","sex", "is.gen54", "is.gen56", "is.batch15", "is.batch22","is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box6", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "comerr8")


# DIFFERENCE BETWEEN FINAL AND INITIAL PREFERENCE (D8 - D1, total and six 5-min bins)
traitcovs[["cpp.diff"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box6", "is.cpp.box7", "is.cpp.box11", "is.cpp.box12", "is.batch2", "is.batch3","is.batch9", "is.batch14", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch21", "is.batch22")
traitcovs[["cpp.diff.p"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box6", "is.cpp.box7", "is.cpp.box11", "is.cpp.box12", "is.batch2", "is.batch3","is.batch9", "is.batch14", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch21", "is.batch22")
traitcovs[["cpp.diff1"]]  <- list("one", "sex","is.gen55", "is.gen56","is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box6", "is.cpp.box7", "is.cpp.box11", "is.cpp.box12", "is.batch2", "is.batch3","is.batch9","is.batch14","is.batch16","is.batch17","is.batch18","is.batch19","is.batch21", "is.batch22")
traitcovs[["cpp.diff2"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box5","is.cpp.box6", "is.cpp.box7", "is.cpp.box11", "is.cpp.box12", "is.batch3", "is.batch14", "is.batch22")
traitcovs[["cpp.diff3"]]  <- list("one", "sex", "is.cpp.box10", "is.batch7", "is.batch15")

traitcovs[["cpp.diff4"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box12", "is.batch10")
traitcovs[["cpp.diff5"]]  <- list("one", "sex", "is.cpp.box3", "is.cpp.box4", "is.cpp.box9", "is.cpp.box12")
traitcovs[["cpp.diff6"]]  <- list("one", "sex", "is.cpp.box2", "is.cpp.box4","is.cpp.box6",  "is.cpp.box11", "is.cpp.box12")


# D1 ACTIVITY (SALINE) - RESPONSE TO NOVELTY (total and six 5-min bins)
traitcovs[["act1.t"]]  <- list("one", "sex", "is.batch14", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "comerr1")
traitcovs[["act1.1"]]  <- list("one", "sex", "is.gen51", "is.batch14", "is.batch20", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "comerr1")
traitcovs[["act1.2"]]  <- list("one", "sex", "is.gen51", "is.batch4", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "comerr1")
traitcovs[["act1.3"]]  <- list("one", "sex", "is.gen52","is.batch7", "is.batch10","is.batch14", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box11", "is.cpp.box12", "comerr1")
traitcovs[["act1.4"]]  <- list("one", "sex", "is.batch10","is.batch14", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9","is.cpp.box10", "is.cpp.box11", "is.cpp.box12", "comerr1")
traitcovs[["act1.5"]]  <- list("one", "sex", "is.batch10","is.batch14", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box9", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12", "comerr1")
traitcovs[["act1.6"]]  <- list("one", "sex", "is.gen51", "is.batch2","is.batch14","is.batch21","is.batch22", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6", "is.cpp.box7", "is.cpp.box8", "is.cpp.box11", "is.cpp.box12", "comerr1")


# D2 ACTIVITY - METH (total and six 5-min bins)
traitcovs[["act2.t"]]  <- list("one", "sex", "is.gen51", "is.gen52","is.gen56", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11","is.cpp.box12", "is.batch8","is.batch15", "is.batch21")
traitcovs[["act2.1"]]  <- list("one", "sex", "is.gen51","is.gen56", "is.batch10", "is.batch15", "is.batch21", "is.cpp.box3", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act2.2"]]  <- list("one", "sex", "is.gen51","is.gen56", "is.batch21", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box12")
traitcovs[["act2.3"]]  <- list("one", "sex", "is.gen51","is.gen56", "is.batch10", "is.batch21", "is.cpp.box7", "is.cpp.box8")
traitcovs[["act2.4"]]  <- list("one", "sex", "is.gen51","is.gen56", "is.batch2", "is.batch4","is.batch15", "is.batch21", "is.cpp.box7", "is.cpp.box8")
traitcovs[["act2.5"]]  <- list("one", "sex","is.batch15",  "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act2.6"]]  <- list("one", "sex", "is.gen52","is.gen53","is.batch10","is.batch12", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10","is.cpp.box11", "is.cpp.box12")


# D3 ACTIVITY - SALINE (total and six 5-min bins)
traitcovs[["act3.t"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55", "is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch14", "is.batch15", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch20", "is.batch21", "is.batch22", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act3.1"]]  <- list("one", "sex", "is.gen51", "is.gen56", "is.batch7", "is.batch14", "is.batch21", "is.batch22", "is.cpp.box2", "is.cpp.box7")
traitcovs[["act3.2"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55","is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch14", "is.batch15", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch20", "is.batch21", "is.batch22", "is.cpp.box7", "is.cpp.box10")
traitcovs[["act3.3"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch14", "is.batch15", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch20", "is.batch21", "is.batch22")
traitcovs[["act3.4"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch14", "is.batch15", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch20", "is.batch21", "is.batch22", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act3.5"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55","is.gen56", "is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch14", "is.batch15", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch20", "is.batch21", "is.batch22", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act3.6"]]  <- list("one", "sex", "is.gen51", "is.gen52", "is.gen53", "is.gen54", "is.gen55", "is.gen56","is.batch2", "is.batch3", "is.batch4", "is.batch5", "is.batch6", "is.batch7", "is.batch8", "is.batch9", "is.batch10", "is.batch11", "is.batch12", "is.batch13", "is.batch14", "is.batch15", "is.batch16", "is.batch17", "is.batch18", "is.batch19", "is.batch20", "is.batch21", "is.batch22", "is.cpp.box3", "is.cpp.box4", "is.cpp.box5", "is.cpp.box6","is.cpp.box7", "is.cpp.box8", "is.cpp.box9","is.cpp.box10", "is.cpp.box11", "is.cpp.box12")


# D4 ACTIVITY - METH (total and six 5-min bins)
traitcovs[["act4.t"]]  <- list("one", "sex", "is.cpp.box7", "is.cpp.box8", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act4.1"]]  <- list("one", "sex", "is.gen52","is.gen56", "is.batch7", "is.batch10", "is.batch14", "is.batch17", "is.cpp.box3", "is.cpp.box5", "is.cpp.box7", "is.cpp.box8", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act4.2"]]  <- list("one", "sex", "is.batch12", "is.cpp.box7")
traitcovs[["act4.3"]]  <- list("one", "sex", "is.cpp.box7")
traitcovs[["act4.4"]]  <- list("one", "sex", "is.batch15",  "is.cpp.box7", "is.cpp.box8")
traitcovs[["act4.5"]]  <- list("one", "sex", "is.cpp.box7", "is.cpp.box8","is.cpp.box9","is.cpp.box10", "is.cpp.box11", "is.cpp.box12")
traitcovs[["act4.6"]]  <- list("one", "sex", "is.gen52","is.gen53",  "is.gen54","is.gen55", "is.batch12", "is.batch17", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box5","is.cpp.box6", "is.cpp.box7", "is.cpp.box8","is.cpp.box9","is.cpp.box10", "is.cpp.box11", "is.cpp.box12")


# D5 ACTIVITY - SALINE (total and six 5-min bins)
traitcovs[["act5.t"]]  <- list("one", "sex", "is.gen51", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.cpp.box7", "is.cpp.box8", "is.cpp.box10", "is.cpp.box11","is.batch6", "is.batch7", "is.batch13", "is.batch15", "is.batch16", "is.batch19", "is.batch20", "is.batch21", "is.batch22")
traitcovs[["act5.1"]]  <- list("one", "sex", "is.gen51", "is.gen53","is.gen55", "is.gen56", "is.cpp.box7","is.cpp.box11","is.batch6", "is.batch7", "is.batch13", "is.batch15",  "is.batch20", "is.batch21", "is.batch22")
traitcovs[["act5.2"]]  <- list("one", "sex", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.cpp.box7", "is.batch7", "is.batch13", "is.batch15", "is.batch20", "is.batch21", "is.batch22")
traitcovs[["act5.3"]]  <- list("one", "sex", "is.gen51", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.cpp.box2", "is.cpp.box4", "is.cpp.box6", "is.batch6", "is.batch7", "is.batch13", "is.batch15", "is.batch16","is.batch17", "is.batch19", "is.batch20", "is.batch22")
traitcovs[["act5.4"]]  <- list("one", "sex", "is.gen51", "is.gen53", "is.gen54", "is.gen55", "is.gen56", "is.cpp.box7", "is.cpp.box8", "is.cpp.box11","is.batch7","is.batch9","is.batch10", "is.batch13", "is.batch15", "is.batch16","is.batch18", "is.batch19", "is.batch20", "is.batch21", "is.batch22")
traitcovs[["act5.5"]]  <- list("one", "sex", "is.gen51", "is.gen54", "is.gen55", "is.gen56","is.cpp.box5",  "is.cpp.box7", "is.cpp.box8","is.cpp.box9", "is.cpp.box10", "is.cpp.box11","is.batch7", "is.batch20", "is.batch22")
traitcovs[["act5.6"]]  <- list("one", "sex", "is.cpp.box3", "is.cpp.box4","is.cpp.box5", "is.cpp.box6","is.cpp.box7", "is.cpp.box8","is.cpp.box9", "is.cpp.box10", "is.cpp.box11","is.cpp.box12",  "is.batch7","is.batch12", "is.batch13")


# D8 ACTIVITY - SALINE (total and six 5-min bins)
traitcovs[["act8.t"]]  <- list("one", "sex","is.gen51","is.gen55", "is.gen56", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box5","is.cpp.box6", "is.cpp.box7", "is.cpp.box8","is.cpp.box9","is.cpp.box11", "is.cpp.box12", "is.batch3", "is.batch7")
traitcovs[["act8.1"]]  <- list("one", "sex", "is.gen56", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4", "is.cpp.box7", "is.cpp.box8","is.cpp.box9","is.cpp.box11", "is.batch3", "is.batch7","is.batch12", "is.batch20", "is.batch21")
traitcovs[["act8.2"]]  <- list("one", "sex","is.gen56", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box6", "is.cpp.box7", "is.cpp.box8","is.cpp.box11", "is.cpp.box12", "is.batch2","is.batch3", "is.batch7","is.batch11","is.batch14","is.batch17")
traitcovs[["act8.3"]]  <- list("one", "sex", "is.cpp.box3", "is.cpp.box4","is.cpp.box7", "is.cpp.box8","is.cpp.box11", "is.cpp.box12", "is.batch2","is.batch3", "is.batch7")
traitcovs[["act8.4"]]  <- list("one", "sex","is.gen51", "is.gen54", "is.gen55", "is.gen56", "is.cpp.box3", "is.cpp.box4", "is.cpp.box7", "is.cpp.box8","is.cpp.box9","is.cpp.box11", "is.cpp.box12", "is.batch15")
traitcovs[["act8.5"]]  <- list("one", "sex","is.gen51","is.gen55", "is.gen56", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box5","is.cpp.box6", "is.cpp.box7", "is.cpp.box8","is.cpp.box9","is.cpp.box10","is.cpp.box11", "is.cpp.box12", "is.batch20")
traitcovs[["act8.6"]]  <- list("one", "sex","is.gen51","is.gen55", "is.gen56", "is.cpp.box2", "is.cpp.box3", "is.cpp.box4","is.cpp.box5","is.cpp.box6", "is.cpp.box7", "is.cpp.box8","is.cpp.box9", "is.cpp.box10","is.cpp.box11", "is.cpp.box12", "is.batch3", "is.batch21")


# LOCOMOTOR SENSITIZATION TO METH (D4 - D2 activity, total and six 5-min bins)
traitcovs[["sens"]]    <- list("one", "sex", "is.gen52", "is.cpp.box7","is.cpp.box10", "is.batch8")
traitcovs[["sens1"]]    <- list("one", "sex", "is.batch4", "is.batch7", "is.batch21")
traitcovs[["sens2"]]    <- list("one", "sex", "is.cpp.box5","is.cpp.box7","is.cpp.box10", "is.cpp.box12")
traitcovs[["sens3"]]    <- list("one", "sex", "is.cpp.box7", "is.batch2", "is.batch11")
traitcovs[["sens4"]]    <- list("one", "sex", "is.gen55","is.gen56", "is.cpp.box7", "is.batch17", "is.batch21")
traitcovs[["sens5"]]    <- list("one", "sex", "is.cpp.box7", "is.batch9", "is.batch11")
traitcovs[["sens6"]]    <- list("one", "sex", "is.cpp.box3","is.cpp.box4","is.cpp.box6","is.cpp.box7","is.cpp.box11","is.cpp.box12")


# OTHER TRAITS
traitcovs[["wild.binary"]]    <- list("one", "sex", "cpp.age")
traitcovs[["tail"]] <-list("one", "sex", "is.gen52", "is.gen53", "is.gen51", "is.gen56", "rip.weight")
traitcovs[["glucose"]] <- list("one", "sex", "glu.weight", "glu.age")
traitcovs[["is.coatA"]] <- list("one", "sex")
traitcovs[["is.coatB"]] <- list("one", "sex")
traitcovs[["is.coatW"]] <- list("one", "sex")



### MAKE SCRIPT TO RUN GEMMA FOR CHOSEN TRAIT ----------------------------------------------

cmds <- c()
for (trait in names(traitcovs)) {

    index.pheno     <- which(pheno.names == trait)
    chosen.covars   <- covars[, unlist(traitcovs[[trait]])]

    write.table(chosen.covars,
                file=paste0("/group/palmer-lab/AIL/qtlmapping/covariates/", trait, ".emp.covs"),
                sep="\t", quote=F, row.names=F, col.names=F)

    for (chrom in 1:19) {
        cmds <- c(cmds, paste0("gemma -g /group/palmer-lab/AIL/GBS/dosage/chr", chrom,
                               ".filtered.dosage.emp -p /group/palmer-lab/AIL/LgSm-DataProcessing/phenos.allgeno.txt -k /group/palmer-lab/AIL/qtlmapping/kinship/chrNot",
                               chrom,".cXX.txt -a /group/palmer-lab/AIL/GBS/dosage/chr",
                               chrom, ".filtered.snpinfo.emp -c /group/palmer-lab/AIL/qtlmapping/covariates/",
                               trait, ".emp.covs -lmm 2 -maf ", MAF, " -o ", trait, ".chr",
                               chrom, " -n ", index.pheno))
    }
}
write.table(cmds, file=paste0("/group/palmer-lab/AIL/code/gemma.alltraits.emp.cmds"),
            row.names=F, col.names=F, quote=F)
