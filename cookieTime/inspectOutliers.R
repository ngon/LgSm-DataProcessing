##### Manual removal of phenotypic outliers
##### 18 Mar 2015

setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/figures/residuals/qq.resid/")
library("ggplot2")
source("../cookieTime/read.pheno.R")
source("../cookieTime/check.normal.quantiles.R")
source("../cookieTime/gg_qq.R")

#### I. LOAD DATA --------------------------------------------------------------

# read phenotype data
alldata <- read.pheno("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles/allData.txt")

# read ppi data with all bins included
allppi <- read.table("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles/ppi.data.ALL.txt", sep="\t", header=T)

info <- read.table("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles/ppi.info.all.txt", sep="\t", header=T)
info <- info[info$id %in% allppi$ID,]

allppi <- data.frame(allppi, info$ppi.weight)

#### II. PLOT RESIDUALS AND IDENTIFY OUTLIERS 

#### A. PREPARE CPP, GLUC, WILD & TAIL PHENOS FOR PLOTTING -------------------
####   All of these phenos have four covariates) 
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
  act8.t = list(pheno='act8.t', cov=c("sex", "gen", "batch", "cpp.box")),

  glucose = list(pheno='glucose', cov=c("sex","gen", "glu.weight", "glu.age")),
  tail = list(pheno="tail", cov=c("sex", "gen", "batch", "rip.weight")))

phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno <- alldata

#### B. IDENTIFY OUTLIERS AND PLOT CPP, GLUC, WILD & TAIL DATA ----------------

# outlierList is a list of lists for each trait. each trait list 
# contains individuals who fall outside the confidence interval specified
# in the call to gg_qq. Their indexes in model$resid are listed in $label.
# $ord.x gives the Observed Values, $z gives the Expected Values, $upper
# and $lower are the bounds for the pointwise confidence interval. 
outlierList=list()

# normQuantiles is a list of lists for each trait that contains the Expected
# and Observed cumulative distributions 
normQuantiles=list()

# THIS LOOP WORKS FOR TRAITS WITH 4 COVARIATES
for (panel in names(panels)) {
  r <- panels[[panel]]
  phenotype <- r$pheno
  covariate <- r$cov
  data        <- pheno[c(phenotype,covariate)]
  names(data) <- c("y","x1", "x2", "x3", "x4")
  model <- lm(y ~ x1 + x2 + x3 + x4, data)
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
  }
# bind pheno ids to new data frame  
pheno.rm.out <- data.frame(pheno$id, pheno.rm.out)


# for each trait, get the ids of animals whose measurements were omitted

out <- list()
for (i in seq_along(outliers)) { out[[i]] <- as.integer(outliers[[i]]) }
names(out) <- names(outliers)
out <- lapply(out, function(x) x[!is.na(x)])

ids <- list()
for (i in seq_along(outliers)){
  m <- subset(pheno[1], row.names(pheno) %in% outliers[[i]])
  ids[[i]] <- m
  }
names(ids) <- names(outliers)

for (i in names(ids)){
lapply(ids[i], write.table, file=paste0(i, "outlierList.txt"), 
       sep="\t", row.names=T, quote=F)
}



#### C. PREPARE PPI PHENOS FOR PLOTTING ------------------------------------

panels <- list(   
  ppi3.logit = list(pheno='ppi3.logit', cov=c("sex","gen", "ppi.box", 
                                       "batch", "ppi.weight")),
  ppi6.logit = list(pheno='ppi6.logit', cov=c("sex", "gen","ppi.box", 
                                       "batch","ppi.weight")),
  ppi12.logit = list(pheno='ppi12.logit', cov=c("sex", "gen","ppi.box", 
                                         "batch","ppi.weight")),
  startle = list(pheno='startle', cov=c("sex", "gen","ppi.box", 
                                      "batch","ppi.weight")),
  habituation = list(pheno='habituation', cov=c("sex","gen", "ppi.box", 
                                          "batch","ppi.weight")))

phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno <- alldata

#### D. IDENTIFY & REMOVE PROBLEMATIC SAMPLES IN PPI DATA -----------------------

# ppi phenotypes do not follow a normal distribution. ppi3, 6 and 12 are
# proportions QTLs will be mapped on transformed values (logit). 
# startle should be transformed with log(10) and habituation is ~normal.

# find the average startle amplitude for the no stimulus trials
allppi$avg.nostim <- (allppi$nostim.1 + allppi$nostim.2 + allppi$nostim.3 + 
                       allppi$nostim.4 + allppi$nostim.5 + allppi$nostim.6 +
                       allppi$nostim.7 + allppi$nostim.8)/8

# find the mice with the smallest differences between avg.nostim and startle.
# i chose 15 as my cutoff for a "suspiciously small difference" and examined
# the data manually. 
allppi$diff.start.nostim <- allppi$startle - allppi$avg.nostim
quantile(allppi$diff.start.nostim, probs=seq(0, 1, 0.05), na.rm=T)
errors <- allppi[allppi$diff.start.nostim <= 15,]

# of the mice whose diff.start.nostim <=15, i'm throwing out mice whose 
# avg.nostim >= 30. many of these appear to be technical errors. i also looked
# at the mice whose diff.start.nostim > 15 and whose avg.nostim >=30, but in
# all of these cases, diff.start.nostim was large enough that an avg.nostim 
# score >=30 did not appear out of the ordinary. notably, 53/59 of the mice w/
# avg.nostim >= 30 were tested in box 3, and a few batches were overrepresented
# in this subset. since i'm looking at raw data, NOT residuals, i don't want to
# be overzealous in removing samples since i expect that covariates will help to
# further control deviation from the mean. 

# i'm also discarding mice for whom the difference between average startle
# and average nostim is <=5.0. i looked at values that were slightly greater
# than 5.0 to see if i'd missed anything, but it appears to be ok. 

# a list of these ids (n=25) are in the file ppi.outliers.txt.
ppiOutliers <- read.table("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles/ppi.outliers.txt", sep="\t", header=T)

data <- pheno[c(1, 129:138)]
matches <- which(data$id %in% ppiOutliers$ppi.outliers)
data.noid <- data[-1]
rm.missing <- lapply(data.noid, function(x) replace(x, matches, "NA"))
rm.missing <- as.numeric(rm.missing)
ppi <- do.call(what=cbind.data.frame, args=rm.missing)
ppi <- data.frame(pheno$id, ppi)
names(ppi)[1] <- "id"

write.table(ppi, "ppi.rm.out.txt", sep="\t", row.names=F, col.names=T)
ppi <-read.table("ppi.rm.out.txt", sep="\t", header=T)

#### E. UPDATE FILES TO BE USED FOR QTL MAPPING -------------------------------
# this is the updated form of Natalia's file 'allData' with covariates and
# phenotypes that have had outliers removed
otherVars <- pheno[c(1:36, )]
pheno.noOutliers <- data.frame(pheno[c(1:36)], pheno.rm.out[-1], 
                           pheno[c(110:116,118:119, 121:122, 124:128)],  
                           ppi.rm.out[-1], pheno[c(152:204)])
write.table(pheno.noOutliers, "C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles/allData.rmOut.txt", col.names=T, row.names=F)

# these are the updated phenotype files for Gemma; they do not include 
# covariates and the header is stored separately.
phenonoOutliers <- data.frame(pheno.noOutliers[c(1, 25:122, 131, 133:134,
                                                  136, 139:149, 157:159)])
write.table(phenonoOutliers, "C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles/phenotypes.orm.txt", sep="\t", col.names=T, row.names=F, quote = F)

