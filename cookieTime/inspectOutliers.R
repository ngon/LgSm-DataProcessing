##### Manual removal of phenotypic outliers
##### 18 Mar 2015

setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles")
library("ggplot2")
source("../cookieTime/read.pheno.R")
source("../cookieTime/binary.from.categorical.R")
source("../cookieTime/check.normal.quantiles.R")
source("../cookieTime/gg_qq.R")


#### I. LOAD DATA AND PREPARE PHENOTYPES -----------------------------------

# read and phenotype data
pheno <- read.pheno("allData.txt")

# make binary phenotypes from some factor variables
binaryPhenos <- cbind(binary.from.categorical(pheno$gen, 
                                              col.names=paste0("is.gen", 50:56)),
                      binary.from.categorical(pheno$coat, 
                                              col.names=paste0("is.coat",
                                                               c("A","B","W"))),
                      binary.from.categorical(pheno$batch, 
                                              col.names=paste0("is.batch", 1:22)),
                      binary.from.categorical(pheno$cpp.box, 
                                              col.names=paste0("is.cpp.box",1:12)),
                      binary.from.categorical(pheno$ppi.box, 
                                              col.names=paste0("is.ppi.box",1:5)))
# combine with phenotypes data frame
pheno <- cbind(pheno, binaryPhenos)

#### II. PLOT RESIDUALS AND IDENTIFY OUTLIERS --------------------------------

### A. PREPARE PHENOS FOR PLOTTING 

### 1. CPP PHENOTYPES (all of which have four covariates) 
panels <- list(
  cpp.diff1 = list(pheno='cpp.diff1', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp.diff2 = list(pheno='cpp.diff2', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp.diff3 = list(pheno='cpp.diff3', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp.diff4 = list(pheno='cpp.diff4', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp.diff5 = list(pheno='cpp.diff5', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp.diff6 = list(pheno='cpp.diff6', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp.diff = list(pheno='cpp.diff', cov=c("sex", "gen", "batch", "cpp.box")),
  
  cpp1.1 = list(pheno='cpp1.1', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp1.2 = list(pheno='cpp1.2', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp1.3 = list(pheno='cpp1.3', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp1.4 = list(pheno='cpp1.4', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp1.5 = list(pheno='cpp1.5', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp1.6 = list(pheno='cpp1.6', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp1.t = list(pheno='cpp1.t', cov=c("sex", "gen", "batch", "cpp.box")),
    
  cpp8.1 = list(pheno='cpp8.1', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp8.2 = list(pheno='cpp8.2', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp8.3 = list(pheno='cpp8.3', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp8.4 = list(pheno='cpp8.4', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp8.5 = list(pheno='cpp8.5', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp8.6 = list(pheno='cpp8.6', cov=c("sex", "gen", "batch", "cpp.box")),
  cpp8.t = list(pheno='cpp8.t', cov=c("sex", "gen", "batch", "cpp.box")),
    
  sc1.1 = list(pheno='sc1.1', cov=c("sex", "gen", "batch", "cpp.box")),
  sc1.2 = list(pheno='sc1.2', cov=c("sex", "gen", "batch", "cpp.box")),
  sc1.3 = list(pheno='sc1.3', cov=c("sex", "gen", "batch", "cpp.box")),
  sc1.4 = list(pheno='sc1.4', cov=c("sex", "gen", "batch", "cpp.box")),
  sc1.5 = list(pheno='sc1.5', cov=c("sex", "gen", "batch", "cpp.box")),
  sc1.6 = list(pheno='sc1.6', cov=c("sex", "gen", "batch", "cpp.box")),
  sc1.t = list(pheno='sc1.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  sc8.1 = list(pheno='sc8.1', cov=c("sex", "gen", "batch", "cpp.box")),
  sc8.2 = list(pheno='sc8.2', cov=c("sex", "gen", "batch", "cpp.box")),
  sc8.3 = list(pheno='sc8.3', cov=c("sex", "gen", "batch", "cpp.box")),
  sc8.4 = list(pheno='sc8.4', cov=c("sex", "gen", "batch", "cpp.box")),
  sc8.5 = list(pheno='sc8.5', cov=c("sex", "gen", "batch", "cpp.box")),
  sc8.6 = list(pheno='sc8.6', cov=c("sex", "gen", "batch", "cpp.box")),
  sc8.t = list(pheno='sc8.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  sens1= list(pheno='sens1', cov=c("sex", "gen", "batch", "cpp.box")),
  sens2= list(pheno='sens2', cov=c("sex", "gen", "batch", "cpp.box")),
  sens3= list(pheno='sens3', cov=c("sex", "gen", "batch", "cpp.box")),
  sens4= list(pheno='sens4', cov=c("sex", "gen", "batch", "cpp.box")),
  sens5= list(pheno='sens5', cov=c("sex", "gen", "batch", "cpp.box")),
  sens6= list(pheno='sens6', cov=c("sex", "gen", "batch", "cpp.box")),
  sens= list(pheno='sens', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act1.1 = list(pheno='act1.1', cov=c("sex","gen", "batch", "cpp.box")),
  act1.2 = list(pheno='act1.2', cov=c("sex", "gen","batch", "cpp.box")),
  act1.3 = list(pheno='act1.3', cov=c("sex", "gen","batch", "cpp.box")),
  act1.4 = list(pheno='act1.4', cov=c("sex","gen", "batch", "cpp.box")),
  act1.5 = list(pheno='act1.5', cov=c("sex","gen", "batch", "cpp.box")),
  act1.6 = list(pheno='act1.6', cov=c("sex","gen", "batch", "cpp.box")),
  act1.t = list(pheno='act1.t', cov=c("sex","gen", "batch", "cpp.box")),
  
  act2.1 = list(pheno='act2.1', cov=c("sex", "gen", "batch", "cpp.box")),
  act2.2 = list(pheno='act2.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act2.3 = list(pheno='act2.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act2.4 = list(pheno='act2.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act2.5 = list(pheno='act2.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act2.6 = list(pheno='act2.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act2.t = list(pheno='act2.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act3.1 = list(pheno='act3.1', cov=c("sex", "gen", "batch", "cpp.box")),
  act3.2 = list(pheno='act3.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act3.3 = list(pheno='act3.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act3.4 = list(pheno='act3.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act3.5 = list(pheno='act3.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act3.6 = list(pheno='act3.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act3.t = list(pheno='act3.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act4.1 = list(pheno='act4.1', cov=c("sex","gen", "batch",  "cpp.box")),
  act4.2 = list(pheno='act4.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act4.3 = list(pheno='act4.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act4.4 = list(pheno='act4.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act4.5 = list(pheno='act4.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act4.6 = list(pheno='act4.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act4.t = list(pheno='act4.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act5.1 = list(pheno='act5.1', cov=c("sex", "gen", "batch", "cpp.box")),
  act5.2 = list(pheno='act5.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act5.3 = list(pheno='act5.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act5.4 = list(pheno='act5.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act5.5 = list(pheno='act5.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act5.6 = list(pheno='act5.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act5.t = list(pheno='act5.t', cov=c("sex", "gen", "batch", "cpp.box")),
  
  act8.1 = list(pheno='act8.1', cov=c("sex", "gen", "batch", "cpp.box")),
  act8.2 = list(pheno='act8.2', cov=c("sex", "gen", "batch", "cpp.box")),
  act8.3 = list(pheno='act8.3', cov=c("sex", "gen", "batch", "cpp.box")),
  act8.4 = list(pheno='act8.4', cov=c("sex", "gen", "batch", "cpp.box")),
  act8.5 = list(pheno='act8.5', cov=c("sex", "gen", "batch", "cpp.box")),
  act8.6 = list(pheno='act8.6', cov=c("sex", "gen", "batch", "cpp.box")),
  act8.t = list(pheno='act8.t', cov=c("sex", "gen", "batch", "cpp.box")))

phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno <- alldata


## Get residual plots and outlier/quantile stats --------------------------

# outlierList is a list of lists for each trait. each trait list 
# contains individuals who fall outside the confidence interval specified
# in the call to gg_qq. Their indexes in model$resid are listed in $label.
# $ord.x gives the Observed Values, $z gives the Expected Values, $upper
# and $lower are the bounds for the pointwise confidence interval. 
outlierList=list()

# normQuantiles is a list of lists for each trait that contains the Expected
# and Observed cumulative distributions 
normQuantiles=list()

# THIS CODE ONLY WORKS FOR TRAITS WITH 4 COVARIATES
for (panel in names(panels)) {
  r <- panels[[panel]]
  phenotype <- r$pheno
  covariate <- r$cov
  data        <- pheno[c(phenotype,covariate)]
  names(data) <- c("y","x1", "x2", "x3", "x4")
  model <- lm(y ~ x1 + x2 +x3 + x4, data)
  x <- (model$resid - mean(model$resid))/sqrt(var(model$resid))
  outlierList[[panel]] <- gg_qq(x, trait=phenotype)
  normQuantiles[[panel]] <- (check.normal.quantiles(model$resid))
}

# identify outliers - they must be both >3 sd from the mean AND fall
# outside the 99% confidence interval. 
outliers=list()
for (panel in names(outlierList)) {
  p <- outlierList[[panel]]  
  x <- as.character(p$ci.out)
  y <- as.character(p$sd3.out)
  outliers[[panel]] <- intersect(x,y)
}

pheno.rm.out <- list()
# remove outlier values and replace with NA
for (panel in names(panels)) {
  r <- panels[[panel]]
  phenotype <- r$pheno
  data <- pheno[phenotype]
  o <- as.integer(outliers[[panel]])
  o <- o[which(complete.cases(o))]
  matches <- which(row.names(data) %in% o)
  rm.missing <- unlist(lapply(data, function(x) replace(x, matches, "NA")))
  rm.missing <- as.numeric(rm.missing)
  pheno.rm.out[[panel]] <- rm.missing
  pheno.rm.out <- do.call(what=cbind.data.frame, args=pheno.rm.out)
  pheno.rm.out <- data.frame(pheno$id, pheno.rm.out)
  }
  





########### ppi traits #################
# gluc = list(pheno='glucose', cov=c("sex", "gluc.weight", "gluc.age")),
# wildness = list(pheno="wild", cov=c("sex", "cpp.age")),
# 
#  panels <- list(     
# pp3 = list(pheno='ppi3.logit', cov=c("sex","gen", "ppi.box", "batch", "ppi.weight")),
# pp6 = list(pheno='ppi6.logit', cov=c("sex", "gen","ppi.box", "batch","ppi.weight")),
# pp12 = list(pheno='ppi12.logit', cov=c("sex", "gen","ppi.box", "batch","ppi.weight")),
# start = list(pheno='startle', cov=c("sex", "gen","ppi.box", "batch","ppi.weight")),
# habit = list(pheno='habituation', cov=c("sex","gen", "ppi.box", "batch","ppi.weight")))


  
 