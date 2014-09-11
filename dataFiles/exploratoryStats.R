source("./oldScripts/SummarySE.R")
library(plyr)
library(ggplot2)


############ ACTIVITY AND CPP DATA #######################

# read in correct data frame
load(cppData, file="./CPPandPPI.11sep14.RData")
data <- cppData

# create activity df for exploratory analysis
# covs are on the left (up to col 14) and traits from 37-119
actWcovs <- cbind(data[c(1:2, 7:10, 13,14)], data[c(37:109,119)])

# paired t tests for activity on meth and saline
senst <- t.test(actWcovs$act2.t, actWcovs$act4.t, paired=TRUE)
senst$p.value
salLoct <- t.test(actWcovs$act3.t, actWcovs$act5.t, paired=TRUE)
salLoct$p.value

# paired t tests for cpp
cpptest <- t.test(actWcovs$cpp8.t, actWcovs$cpp1.t, paired=TRUE)
cpptest$p.value

# create cpp.diff.sec
actWcovs$cpp.diff.sec <- actWcovs$cpp8.t - actWcovs$cpp1.t

#################### MALES AND FEMALES ################################

# create male and female  dfs
actF <- actWcovs[actWcovs$sex == "F",]
actM <- actWcovs[actWcovs$sex == "M",]

# t tests on male and female activity (all are p 0.02(sz) or lower)
d1at <- t.test(actF$act1.t, actM$act1.t, paired=F)
d1at$p.value
d2at <- t.test(actF$act2.t, actM$act2.t, paired=F)
d2at$p.value
d3at <- t.test(actF$act3.t, actM$act3.t, paired=F)
d3at$p.value
d4at <- t.test(actF$act4.t, actM$act4.t, paired=F)
d4at$p.value
d5at <- t.test(actF$act5.t, actM$act5.t, paired=F)
d5at$p.value
d8at <- t.test(actF$act8.t, actM$act8.t, paired=F)
d8at$p.value
senMFt <- t.test(actF$sens, actM$sens, paired=F)
senMFt$p.value
CPPmft <- t.test(actF$cpp8.t, actM$cpp1.t, paired=F)
CPPmft$p.value
CPPmft5min <- t.test(actF$cpp8.1, actM$cpp1.1, paired=F)
CPPmft5min$p.value

# not significant
cppDiff <- t.test(actF$cpp.diff.sec, actM$cpp.diff.sec, paired=F)
cppDiff$p.value # 0.227



########### FACTOR COVARIATE ANALYSIS ###################
act <- actWcovs 

sensGen <- lm(act$sens ~ as.factor(act$gen)); summary(sensGen) # gen52, p=0.0303
sensSex <- lm(act$sens ~ act$sex); summary(sensSex) # sex, p=0.0167
sensBox <- lm(act$sens ~ as.factor(act$cpp.box)); summary(sensBox) # box7, p=0.006
sensBatch <- lm(act$sens ~ as.factor(act$batch)); summary(sensBatch) # batch 8

d1Gen <- lm(act$act1.t ~ as.factor(act$gen)); summary(d1Gen) # no effect
d1Box <- lm(act$act1.t ~ as.factor(act$cpp.box)); summary(d1Box) # effects of box 3-9,11,12
d1Batch <- lm(act$act1.t ~ as.factor(act$batch)); summary(d1Batch) # batch 14

d2Box <- lm(act$act2.t ~ as.factor(act$cpp.box)); summary(d2Box) # box 5,7,8,10,12
d2Gen <- lm(act$act2.t ~ as.factor(act$gen)); summary(d2Gen) # gen 52, 56
d2Batch <- lm(act$act2.t ~ as.factor(act$batch)); summary(d2Batch) # 8

d3Box <- lm(act$act3.t ~ as.factor(act$cpp.box)); summary(d3Box) # box 7,8,10,11,12
d3Gen <- lm(act$act3.t ~ as.factor(act$gen)); summary(d3Gen) # all gen except 52
d3Batch <- lm(act$act3.t ~ as.factor(act$batch)); summary(d3Batch) # all batches

d4Box <- lm(act$act4.t ~ as.factor(act$cpp.box)); summary(d4Box) # box 7,8,11
d4Gen <- lm(act$act4.t ~ as.factor(act$gen)); summary(d4Gen) # no effect
d4Batch <- lm(act$act4.t ~ as.factor(act$batch)); summary(d4Batch) # no effect

d5Box <- lm(act$act5.t ~ as.factor(act$cpp.box)); summary(d5Box) # box 7,8
d5Gen <- lm(act$act5.t ~ as.factor(act$gen)); summary(d5Gen) # gen 51, 53-56
d5Batch <- lm(act$act5.t ~ as.factor(act$batch)); summary(d5Batch) # batch 7,13,15,20-22

d8Box <- lm(act$act8.t ~ as.factor(act$cpp.box)); summary(d8Box) # box 3,4,6-9,11,12
d8Gen <- lm(act$act8.t ~ as.factor(act$gen)); summary(d8Gen) # gen 56
d8Batch <- lm(act$act8.t ~ as.factor(act$batch)); summary(d8Batch) # 3,7

cppBox <- lm(act$cpp8.t ~ as.factor(act$cpp.box)); summary(cppBox) # all except box 5
cppGen <- lm(act$cpp8.t ~ as.factor(act$gen)); summary(cppGen) # gen 53,56
cppBatch <- lm(act$cpp8.t ~ as.factor(act$batch)); summary(cppBatch) # 2,10,13,15,21

ippBox <- lm(act$cpp1.t ~ as.factor(act$cpp.box)); summary(ippBox) # box 3,8-10
ippGen <- lm(act$cpp1.t ~ as.factor(act$gen)); summary(ippGen) # gen 54,56
ippBatch <- lm(act$cpp1.t ~ as.factor(act$batch)); summary(ippBatch) # 13,14,17,21,22


################### QUANTATIVE COVARIATES #######################

## no significant correlations between quantitative covariates and any of the traits
## below that are not already encompassed by another covariate (e.g. sex&weight)

panels <-
        list(               
                ipp.ratio= list (pheno="cpp1.t",   cov="mf.ratio",    col="orange"),
                cpp.ratio= list (pheno="cpp8.t",   cov="mf.ratio",    col="orange"),
                cppdiff.ratio= list (pheno="cpp.diff.sec", cov="mf.ratio",    col="orange"),
                act1.t.ratio= list (pheno="act1.t",   cov="mf.ratio",    col="orange"),
                act2.t.ratio= list (pheno="act2.t",   cov="mf.ratio",    col="orange"),
                act3.t.ratio= list (pheno="act3.t",   cov="mf.ratio",    col="orange"),
                act4.t.ratio= list (pheno="act4.t",   cov="mf.ratio",    col="orange"),
                act5.t.ratio= list (pheno="act5.t",   cov="mf.ratio",    col="orange"),
                act8.t.ratio= list (pheno="act8.t",   cov="mf.ratio",    col="orange"),
                sens.ratio= list (pheno="sens",   cov="mf.ratio",    col="orange"),
                
                ipp.weight= list (pheno="cpp1.t",   cov="d8.weight",    col="green"),
                cpp.weight= list (pheno="cpp8.t",   cov="d8.weight",    col="green"),
                cppdiff.weight= list (pheno="cpp.diff.sec", cov="d8.weight",    col="green"),
                
                
                ipp.sc= list (pheno="cpp1.t",   cov="sc1.t",    col="dodgerblue"),
                cpp.sc= list (pheno="cpp8.t",   cov="sc8.t",    col="dodgerblue"),
                cppdiff.8sc= list (pheno="cpp.diff.sec", cov="sc8.t",    col="dodgerblue"),
                cppdiff.1sc= list (pheno="cpp.diff.sec", cov="sc1.t",    col="dodgerblue"),
                act1.t.sc= list (pheno="act1.t",   cov="sc1.t",    col="dodgerblue"),
                act8.t.sc= list (pheno="act8.t",   cov="sc8.t",    col="dodgerblue")
        )



phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno      <- act
n     <- nrow(pheno)
pheno <- transform(pheno,mf.ratio = mf.ratio + rnorm(n,sd = 0.5))



for (panel in names(panels)) {      
        r         <- panels[[panel]]
        phenotype <- r$pheno
        covariate <- r$cov
        panel.col <- r$col
        data        <- pheno[c(phenotype,covariate)]
        #data <- cbind(data, pheno)
        names(data) <- c("y","x")
        model <- lm(y ~ x,data)
        pve   <- summary(model)$r.squared
        
        filename=paste0(covariate,"_", phenotype, ".png", sep="")
        png(file=filename, height=300,width=300)
        
        print(ggplot(data, aes(x=x,y=y)) +
                      geom_point(color=panel.col)+
                      geom_smooth(method=lm, se=TRUE) +
                      xlab(paste0(covariate, " (PVE= ", round(100*pve, digits=2), "%)"))+
                      ylab(paste0(phenotype))+
                      guides(color=FALSE)+
                      theme_bw()+
                      theme(axis.title.x = element_text(colour="black", size=14),
                            axis.text = element_text(colour="black", size=12),
                            axis.title.y = element_text(colour="black", size=14),
                            plot.title = element_text(colour="black", size=12),
                            legend.position="none"))
        
        dev.off()      
}












