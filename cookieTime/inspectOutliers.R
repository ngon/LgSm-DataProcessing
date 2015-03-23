##### Manual removal of phenotypic outliers
##### 18 Mar 2015

setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles")
library("ggplot2")
source("../cookieTime/read.pheno.R")
source("../cookieTime/binary.from.categorical.R")
source("../cookieTime/gg_qq.R")

#### CONTENTS
# 
### 1. LOAD DATA AND PREPARE PHENOTYPES  
### II. PLOT RESIDUALS AND IDENTIFY OUTLIERS
###         A. CPP/LOCOMOTOR PHENOTYPES
###         C. PPI PHENOTYPES
###         D. OTHER PHENOTYPES
### 3. REMOVE OUTLIERS

#####################################################################################
### 1. LOAD DATA AND PREPARE PHENOTYPES ###

# read and phenotype data
pheno <- read.pheno("allData.txt")

# make binary phenotypes from some factor variables
binaryPhenos <- cbind(binary.from.categorical(pheno$gen, 
                                              col.names=paste0("is.gen", 50:56)),
                      binary.from.categorical(pheno$coat, 
                                              col.names=paste0("is.coat", c("A","B","W"))),
                      binary.from.categorical(pheno$batch, 
                                              col.names=paste0("is.batch", 1:22)),
                      binary.from.categorical(pheno$cpp.box, 
                                              col.names=paste0("is.cpp.box",1:12)),
                      binary.from.categorical(pheno$ppi.box, 
                                              col.names=paste0("is.ppi.box",1:5)))

# combine with phenotypes data frame
pheno <- cbind(pheno, binaryPhenos)


#####################################################################################
### II. PLOT RESIDUALS AND IDENTIFY OUTLIERS ###

### A. PREPARE PHENOS FOR PLOTTING 

### 1. CPP PHENOTYPES (all of which have four covariates) 
panels <- list(
  cppDiff1 = list(pheno='cpp.diff1', cov=c("sex", "gen", "batch", "cpp.box")),
  cppDiff2 = list(pheno='cpp.diff2', cov=c("sex", "gen", "batch", "cpp.box")),
  cppDiff3 = list(pheno='cpp.diff3', cov=c("sex", "gen", "batch", "cpp.box")),
  cppDiff4 = list(pheno='cpp.diff4', cov=c("sex", "gen", "batch", "cpp.box")),
  cppDiff5 = list(pheno='cpp.diff5', cov=c("sex", "gen", "batch", "cpp.box")),
  cppDiff6 = list(pheno='cpp.diff6', cov=c("sex", "gen", "batch", "cpp.box")),
  cppDiff = list(pheno='cpp.diff', cov=c("sex", "gen", "batch", "cpp.box")),
  
  ipp11 = list(pheno='cpp1.1', cov=c("sex", "gen", "batch", "cpp.box")),
  ipp12 = list(pheno='cpp1.2', cov=c("sex", "gen", "batch", "cpp.box")),
  ipp13 = list(pheno='cpp1.3', cov=c("sex", "gen", "batch", "cpp.box")),
  ipp14 = list(pheno='cpp1.4', cov=c("sex", "gen", "batch", "cpp.box")),
  ipp15 = list(pheno='cpp1.5', cov=c("sex", "gen", "batch", "cpp.box")),
  ipp16 = list(pheno='cpp1.6', cov=c("sex", "gen", "batch", "cpp.box")),
  ipp1T = list(pheno='cpp1.t', cov=c("sex", "gen", "batch", "cpp.box")),
    
  cpp81 = list(pheno='cpp8.1', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp82 = list(pheno='cpp8.2', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp83 = list(pheno='cpp8.3', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp84 = list(pheno='cpp8.4', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp85 = list(pheno='cpp8.5', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp86 = list(pheno='cpp8.6', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp8T = list(pheno='cpp8.t', cov=c("sex", "gen", "batch", "cpp.box")),
    
  sc11 = list(pheno='sc1.1', cov=c("sex", "gen", "batch", "cpp.box")),
  sc12 = list(pheno='sc1.2', cov=c("sex", "gen", "batch", "cpp.box")),
  sc13 = list(pheno='sc1.3', cov=c("sex", "gen", "batch", "cpp.box")),
  sc14 = list(pheno='sc1.4', cov=c("sex", "gen", "batch", "cpp.box")),
  sc15 = list(pheno='sc1.5', cov=c("sex", "gen", "batch", "cpp.box")),
  sc16 = list(pheno='sc1.6', cov=c("sex", "gen", "batch", "cpp.box")),
  sc1T = list(pheno='sc1.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  sc81 = list(pheno='sc8.1', cov=c("sex", "gen", "batch", "cpp.box")),
  sc82 = list(pheno='sc8.2', cov=c("sex", "gen", "batch", "cpp.box")),
  sc83 = list(pheno='sc8.3', cov=c("sex", "gen", "batch", "cpp.box")),
  sc84 = list(pheno='sc8.4', cov=c("sex", "gen", "batch", "cpp.box")),
  sc85 = list(pheno='sc8.5', cov=c("sex", "gen", "batch", "cpp.box")),
  sc86 = list(pheno='sc8.6', cov=c("sex", "gen", "batch", "cpp.box")),
  sc8T = list(pheno='sc8.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  sens1= list(pheno='sens1', cov=c("sex", "gen", "batch", "cpp.box")),
  sens2= list(pheno='sens2', cov=c("sex", "gen", "batch", "cpp.box")),
  sens3= list(pheno='sens3', cov=c("sex", "gen", "batch", "cpp.box")),
  sens4= list(pheno='sens4', cov=c("sex", "gen", "batch", "cpp.box")),
  sens5= list(pheno='sens5', cov=c("sex", "gen", "batch", "cpp.box")),
  sens6= list(pheno='sens6', cov=c("sex", "gen", "batch", "cpp.box")),
  sens= list(pheno='sens', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act11 = list(pheno='act1.1', cov=c("sex","gen", "batch", "cpp.box")),
  act12 = list(pheno='act1.2', cov=c("sex", "gen","batch", "cpp.box")),
  act13 = list(pheno='act1.3', cov=c("sex", "gen","batch", "cpp.box")),
  act14 = list(pheno='act1.4', cov=c("sex","gen", "batch", "cpp.box")),
  act15 = list(pheno='act1.5', cov=c("sex","gen", "batch", "cpp.box")),
  act16 = list(pheno='act1.6', cov=c("sex","gen", "batch", "cpp.box")),
  act1T = list(pheno='act1.t', cov=c("sex","gen", "batch", "cpp.box")),
  
  act21 = list(pheno='act2.1', cov=c("sex", "gen", "batch", "cpp.box")),
  act22 = list(pheno='act2.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act23 = list(pheno='act2.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act24 = list(pheno='act2.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act25 = list(pheno='act2.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act26 = list(pheno='act2.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act2T = list(pheno='act2.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act31 = list(pheno='act3.1', cov=c("sex", "gen", "batch", "cpp.box")),
  act32 = list(pheno='act3.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act33 = list(pheno='act3.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act34 = list(pheno='act3.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act35 = list(pheno='act3.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act36 = list(pheno='act3.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act3T = list(pheno='act3.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act41 = list(pheno='act4.1', cov=c("sex","gen", "batch",  "cpp.box")),
  act42 = list(pheno='act4.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act43 = list(pheno='act4.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act44 = list(pheno='act4.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act45 = list(pheno='act4.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act46 = list(pheno='act4.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act4T = list(pheno='act4.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act51 = list(pheno='act5.1', cov=c("sex", "gen", "batch", "cpp.box")),
  act52 = list(pheno='act5.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act53 = list(pheno='act5.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act54 = list(pheno='act5.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act55 = list(pheno='act5.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act56 = list(pheno='act5.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act5T = list(pheno='act5.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act81 = list(pheno='act8.1', cov=c("sex", "gen", "batch", "cpp.box")),
  act82 = list(pheno='act8.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act83 = list(pheno='act8.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act84 = list(pheno='act8.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act85 = list(pheno='act8.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act86 = list(pheno='act8.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act8T = list(pheno='act8.t', cov=c("sex", "gen", "batch", "cpp.box")))

 phenotypes <- unique(sapply(panels,function(x)x$pheno))


############### initial place preference ##########################################

  r <- panels[["ipp1T"]]  # change this for each phenotype name
  phenotype <- r$pheno
  covariate <- r$cov
  data        <- pheno[c(phenotype,covariate)]
  names(data) <- c("y","x1", "x2", "x3", "x4")
  model <- lm(y ~ x1 + x2 +x3 + x4, data)
  x <- rstudent(model)
ipp1T.out <- gg_qq(x, trait=phenotype)


r <- panels[["ipp11"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
ipp11.out <-gg_qq(x, trait=phenotype)

r <- panels[["ipp12"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)

ipp12.out <-gg_qq(x, trait=phenotype)

r <- panels[["ipp13"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)

ipp13.out <-gg_qq(x, trait=phenotype)

r <- panels[["ipp14"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)

ipp14.out <-gg_qq(x, trait=phenotype)

r <- panels[["ipp15"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
ipp15.out <-gg_qq(x, trait=phenotype)

r <- panels[["ipp16"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
ipp16.out <-gg_qq(x, trait=phenotype)

############ cpp8 phenos ################################################################

r <- panels[["cpp8T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cpp8T.out <-gg_qq(x, trait=phenotype)

r <- panels[["cpp81"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cpp81.out <-gg_qq(x, trait=phenotype)

r <- panels[["cpp82"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cpp82.out <-gg_qq(x, trait=phenotype)

r <- panels[["cpp83"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cpp83.out <-gg_qq(x, trait=phenotype)

r <- panels[["cpp84"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cpp84.out <-gg_qq(x, trait=phenotype)

r <- panels[["cpp85"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cpp85.out <-gg_qq(x, trait=phenotype)

r <- panels[["cpp86"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cpp86.out <-gg_qq(x, trait=phenotype)


################# cpp diff plots #########################################################

r <- panels[["cppDiff"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cppDiff.out <-gg_qq(x, trait=phenotype)

r <- panels[["cppDiff1"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cppDiff1.out <-gg_qq(x, trait=phenotype)

r <- panels[["cppDiff2"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cppDiff2.out <-gg_qq(x, trait=phenotype)

r <- panels[["cppDiff3"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cppDiff3.out <-gg_qq(x, trait=phenotype)

r <- panels[["cppDiff4"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cppDiff4.out <-gg_qq(x, trait=phenotype)

r <- panels[["cppDiff5"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cppDiff5.out <-gg_qq(x, trait=phenotype)

r <- panels[["cppDiff6"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
cppDiff6.out <-gg_qq(x, trait=phenotype)



################### act day 1 phenotypes #################################################

r <- panels[["act1T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act1T.out <-gg_qq(x, trait=phenotype)

r <- panels[["act11"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act11.out <-gg_qq(x, trait=phenotype)

r <- panels[["act12"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act12.out <-gg_qq(x, trait=phenotype)

r <- panels[["act13"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act13.out <-gg_qq(x, trait=phenotype)

r <- panels[["act14"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act14.out <-gg_qq(x, trait=phenotype)

r <- panels[["act15"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act15.out <-gg_qq(x, trait=phenotype)

r <- panels[["act16"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act16.out <-gg_qq(x, trait=phenotype)

################### act day 2 phenotypes #################################################

r <- panels[["act2T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act2T.out <-gg_qq(x, trait=phenotype)

r <- panels[["act21"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act21.out <-gg_qq(x, trait=phenotype)

r <- panels[["act22"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act22.out <-gg_qq(x, trait=phenotype)

r <- panels[["act23"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act23.out <-gg_qq(x, trait=phenotype)

r <- panels[["act24"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act24.out <-gg_qq(x, trait=phenotype)

r <- panels[["act25"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act25.out <-gg_qq(x, trait=phenotype)

r <- panels[["act26"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act26.out <-gg_qq(x, trait=phenotype)

################### act day 3 phenotypes #################################################

r <- panels[["act3T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act3T.out <-gg_qq(x, trait=phenotype)

r <- panels[["act31"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act31.out <-gg_qq(x, trait=phenotype)

r <- panels[["act32"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act32.out <-gg_qq(x, trait=phenotype)

r <- panels[["act33"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act33.out <-gg_qq(x, trait=phenotype)

r <- panels[["act34"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act34.out <-gg_qq(x, trait=phenotype)

r <- panels[["act35"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act35.out <-gg_qq(x, trait=phenotype)

r <- panels[["act36"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act36.out <-gg_qq(x, trait=phenotype)

################### act day 4 phenotypes #################################################

r <- panels[["act4T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act4T.out <-gg_qq(x, trait=phenotype)

r <- panels[["act41"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act41.out <-gg_qq(x, trait=phenotype)

r <- panels[["act42"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act42.out <-gg_qq(x, trait=phenotype)

r <- panels[["act43"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act43.out <-gg_qq(x, trait=phenotype)

r <- panels[["act44"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act44.out <-gg_qq(x, trait=phenotype)

r <- panels[["act45"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act45.out <-gg_qq(x, trait=phenotype)

r <- panels[["act46"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act46.out <-gg_qq(x, trait=phenotype)

################### act day 5 phenotypes #################################################

r <- panels[["act5T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act5T.out <-gg_qq(x, trait=phenotype)

r <- panels[["act51"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act51.out <-gg_qq(x, trait=phenotype)

r <- panels[["act52"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act52.out <-gg_qq(x, trait=phenotype)

r <- panels[["act53"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act53.out <-gg_qq(x, trait=phenotype)

r <- panels[["act54"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act54.out <-gg_qq(x, trait=phenotype)

r <- panels[["act55"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act55.out <-gg_qq(x, trait=phenotype)

r <- panels[["act56"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act56.out <-gg_qq(x, trait=phenotype)

################### act day 8 phenotypes #################################################

r <- panels[["act8T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act8T.out <-gg_qq(x, trait=phenotype)

r <- panels[["act81"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act81.out <-gg_qq(x, trait=phenotype)

r <- panels[["act82"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act82.out <-gg_qq(x, trait=phenotype)

r <- panels[["act83"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act83.out <-gg_qq(x, trait=phenotype)

r <- panels[["act84"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act84.out <-gg_qq(x, trait=phenotype)

r <- panels[["act85"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act85.out <-gg_qq(x, trait=phenotype)

r <- panels[["act86"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
act86.out <-gg_qq(x, trait=phenotype)


################### side changes day 1 #################################################

r <- panels[["sc1T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc1T.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc11"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc11.out <-gg_qq(x, trait=phenotype)

r <- panels[["sc12"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc12.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc13"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc13.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc14"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc14.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc15"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc15.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc16"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc16.out <-gg_qq(x, trait=phenotype)


################### side changes day 8 #################################################

r <- panels[["sc8T"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc8T.out <-gg_qq(x, trait=phenotype)

r <- panels[["sc81"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc81.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc82"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc82.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc83"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc83.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc84"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc84.out <-gg_qq(x, trait=phenotype)


r <- panels[["sc85"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc85.out <-gg_qq(x, trait=phenotype)

r <- panels[["sc86"]]  # change this for each phenotype name
phenotype <- r$pheno
covariate <- r$cov
data        <- pheno[c(phenotype,covariate)]
names(data) <- c("y","x1", "x2", "x3", "x4")
model <- lm(y ~ x1 + x2 +x3 + x4, data)
x <- rstudent(model)
sc86.out <-gg_qq(x, trait=phenotype)



########### ppi traits #################



# gluc = list(pheno='glucose', cov=c("sex")),
# wildness = list(pheno="wild", cov=c("sex")),
# 
#  panels <- list(     
# pp3 = list(pheno='ppi3.logit', cov=c("sex","gen", "ppi.box", "batch", "ppi.weight")),
# pp6 = list(pheno='ppi6.logit', cov=c("sex", "gen","ppi.box", "batch","ppi.weight")),
# pp12 = list(pheno='ppi12.logit', cov=c("sex", "gen","ppi.box", "batch","ppi.weight")),
# start = list(pheno='startle', cov=c("sex", "gen","ppi.box", "batch","ppi.weight")),
# habit = list(pheno='habituation', cov=c("sex","gen", "ppi.box", "batch","ppi.weight")))


  
 