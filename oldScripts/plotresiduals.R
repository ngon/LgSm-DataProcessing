
## Examining residuals, correlations for each phenotype model ##
library(plyr)

## 1. load data and change a few numbers to factors
alldata <- read.table("./dataFiles//allData.txt", 
                   sep="\t", header=TRUE, na.strings = "NA")


# classes<- lapply(covs, class)
alldata$batch <- as.factor(alldata$batch)
alldata$gen <- as.factor(alldata$gen)
alldata$cpp.box <- as.factor(alldata$cpp.box)
alldata$ppi.box <- as.factor(alldata$ppi.box)


## 2. Linear models

## cpp

panels <- list(
#         cppDiff = list(pheno='cpp.diff', cov=c("sex", "gen", "batch", "cpp.box")),
#         
#         ipp11 = list(pheno='cpp1.1', cov=c("sex", "gen", "batch", "cpp.box")),
#         ipp12 = list(pheno='cpp1.2', cov=c("sex", "gen", "batch", "cpp.box")),
#         ipp13 = list(pheno='cpp1.3', cov=c("sex", "gen", "batch", "cpp.box")),
#         ipp14 = list(pheno='cpp1.4', cov=c("sex", "gen", "batch", "cpp.box")),
#         ipp15 = list(pheno='cpp1.5', cov=c("sex", "gen", "batch", "cpp.box")),
#         ipp16 = list(pheno='cpp1.6', cov=c("sex", "gen", "batch", "cpp.box")),
#         ipp1T = list(pheno='cpp1.t', cov=c("sex", "gen", "batch", "cpp.box")),
#         ipp1P = list(pheno='cpp1.p', cov=c("sex", "gen", "batch", "cpp.box")),
#         
#         cpp81 = list(pheno='cpp8.1', cov=c("sex", "gen", "batch", "cpp.box")),
#         cpp82 = list(pheno='cpp8.2', cov=c("sex", "gen", "batch", "cpp.box")),
#         cpp83 = list(pheno='cpp8.3', cov=c("sex", "gen", "batch", "cpp.box")),
#         cpp84 = list(pheno='cpp8.4', cov=c("sex", "gen", "batch", "cpp.box")),
#         cpp85 = list(pheno='cpp8.5', cov=c("sex", "gen", "batch", "cpp.box")),
#         cpp86 = list(pheno='cpp8.6', cov=c("sex", "gen", "batch", "cpp.box")),
#         cpp8T = list(pheno='cpp8.t', cov=c("sex", "gen", "batch", "cpp.box")),
#         cpp8P = list(pheno='cpp8.p', cov=c("sex", "gen", "batch", "cpp.box")),
#         
#         sc11 = list(pheno='sc1.1', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc12 = list(pheno='sc1.2', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc13 = list(pheno='sc1.3', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc14 = list(pheno='sc1.4', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc15 = list(pheno='sc1.5', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc16 = list(pheno='sc1.6', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc1T = list(pheno='sc1.t', cov=c("sex", "gen", "batch", "cpp.box")),
#         
#         sc81 = list(pheno='sc8.1', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc82 = list(pheno='sc8.2', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc83 = list(pheno='sc8.3', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc84 = list(pheno='sc8.4', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc85 = list(pheno='sc8.5', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc86 = list(pheno='sc8.6', cov=c("sex", "gen", "batch", "cpp.box")),
#         sc8T = list(pheno='sc8.t', cov=c("sex", "gen", "batch", "cpp.box")),
#         
#         sensitization = list(pheno='sens', cov=c("sex", "gen", "batch", "cpp.box")),
        
        
#         gluc = list(pheno='glucose', cov=c("sex")),
#         wildness = list(pheno="wild", cov=c("sex")),
        
#         pp3 = list(pheno='ppi3', cov=c("sex", "ppi.box", "batch")),
#         pp6 = list(pheno='ppi6', cov=c("sex", "ppi.box", "batch")),
#         pp12 = list(pheno='ppi12', cov=c("sex", "ppi.box", "batch")),
#         start = list(pheno='startle', cov=c("sex", "ppi.box", "batch")),
#         habit = list(pheno='habituation', cov=c("sex", "ppi.box", "batch")),
        
#         
#         act11 = list(pheno='act1.1', cov=c("sex", "batch", "cpp.box")),
#         act12 = list(pheno='act1.2', cov=c("sex", "batch", "cpp.box")),
#         act13 = list(pheno='act1.3', cov=c("sex", "batch", "cpp.box")),
#         act14 = list(pheno='act1.4', cov=c("sex", "batch", "cpp.box")),
#         act15 = list(pheno='act1.5', cov=c("sex", "batch", "cpp.box")),
#         act16 = list(pheno='act1.6', cov=c("sex", "batch", "cpp.box")),
#         act1T = list(pheno='act1.t', cov=c("sex", "batch", "cpp.box")),
        
#         act21 = list(pheno='act2.1', cov=c("sex", "gen", "batch", "cpp.box")),
#         act22 = list(pheno='act2.2', cov=c("sex", "gen", "batch", "cpp.box")),
#         act23 = list(pheno='act2.3', cov=c("sex", "gen", "batch", "cpp.box")),
#         act24 = list(pheno='act2.4', cov=c("sex", "gen", "batch", "cpp.box")),
#         act25 = list(pheno='act2.5', cov=c("sex", "gen", "batch", "cpp.box")),
#         act26 = list(pheno='act2.6', cov=c("sex", "gen", "batch", "cpp.box")),
#         act2T = list(pheno='act2.t', cov=c("sex", "gen", "batch", "cpp.box")),
#         
#         act31 = list(pheno='act3.1', cov=c("sex", "gen", "batch", "cpp.box")),
#         act32 = list(pheno='act3.2', cov=c("sex", "gen", "batch", "cpp.box")),
#         act33 = list(pheno='act3.3', cov=c("sex", "gen", "batch", "cpp.box")),
#         act34 = list(pheno='act3.4', cov=c("sex", "gen", "batch", "cpp.box")),
#         act35 = list(pheno='act3.5', cov=c("sex", "gen", "batch", "cpp.box")),
#         act36 = list(pheno='act3.6', cov=c("sex", "gen", "batch", "cpp.box")),
#         act3T = list(pheno='act3.t', cov=c("sex", "gen", "batch", "cpp.box")),
        
        act41 = list(pheno='act4.1', cov=c("sex", "cpp.box")),
        act42 = list(pheno='act4.2', cov=c("sex", "cpp.box")),
        act43 = list(pheno='act4.3', cov=c("sex", "cpp.box")),
        act44 = list(pheno='act4.4', cov=c("sex", "cpp.box")),
        act45 = list(pheno='act4.5', cov=c("sex", "cpp.box")),
        act46 = list(pheno='act4.6', cov=c("sex", "cpp.box")),
        act4T = list(pheno='act4.t', cov=c("sex", "cpp.box")),
        
#         act51 = list(pheno='act5.1', cov=c("sex", "gen", "batch", "cpp.box")),
#         act52 = list(pheno='act5.2', cov=c("sex", "gen", "batch", "cpp.box")),
#         act53 = list(pheno='act5.3', cov=c("sex", "gen", "batch", "cpp.box")),
#         act54 = list(pheno='act5.4', cov=c("sex", "gen", "batch", "cpp.box")),
#         act55 = list(pheno='act5.5', cov=c("sex", "gen", "batch", "cpp.box")),
#         act56 = list(pheno='act5.6', cov=c("sex", "gen", "batch", "cpp.box")),
#         act5T = list(pheno='act5.t', cov=c("sex", "gen", "batch", "cpp.box")),
#         
#         act81 = list(pheno='act8.1', cov=c("sex", "gen", "batch", "cpp.box")),
#         act82 = list(pheno='act8.2', cov=c("sex", "gen", "batch", "cpp.box")),
#         act83 = list(pheno='act8.3', cov=c("sex", "gen", "batch", "cpp.box")),
#         act84 = list(pheno='act8.4', cov=c("sex", "gen", "batch", "cpp.box")),
#         act85 = list(pheno='act8.5', cov=c("sex", "gen", "batch", "cpp.box")),
#         act86 = list(pheno='act8.6', cov=c("sex", "gen", "batch", "cpp.box")),
#         act8T = list(pheno='act8.t', cov=c("sex", "gen", "batch", "cpp.box")))


phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno      <- phenoData


for (panel in names(panels)) { 
        
        # Get the phenotype and covariate to investigate, and the colour for
        # drawing the points in the scatterplot.
        r         <- panels[[panel]]
        phenotype <- r$pheno
        covariate <- r$cov
     
        data        <- pheno[c(phenotype,covariate)]
        names(data) <- c("y","x1", "x2")
        model <- lm(y ~ x1 + x2, data)
        resid <- model$residuals
        
        filename=paste0(covariate,"_", phenotype, ".png", sep="")
        png(file=filename)
        
        plot(resid, ylab=paste0(phenotype))
        
        
#         print(ggplot(data, aes(x,y)) +
#                       geom_point(color=panel.col)+
#                       geom_smooth(method=lm, se=TRUE) +
#                       xlab(paste0(covariate, " (PVE= ", round(100*pve, digits=2), "%"))+
#                       ylab(paste0(phenotype))+
#                       theme(axis.title.x = element_text(colour="black", size=12),
#                             axis.text = element_text(colour="black", size=8),
#                             axis.title.y = element_text(colour="black", size=8),
#                             plot.title = element_text(colour="black", size=8))) 
        
        dev.off()      
}

png(file="./sex_glucose.png")
plot(glu$resid, ylab="glucose")
dev.off()

png(file="./sex_wildness.png")
plot(wildness$resid, ylab="wildness")
dev.off()
