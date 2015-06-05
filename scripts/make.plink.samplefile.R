### MAKE SAMPLE FILE FOR PLINK v1.9 ###
### Sample file is in 'oxford format', the format used in IMPUTE2 and related
### tools. See http://www.stats.ox.ac.uk/~marchini/software/gwas/file_format.html#Sample_File_Format_ for details.

geno.samples<- read.table("./genotyped.samples.txt", as.is=T, header=F)
phenos <- read.table("./dataFiles/phenotypes/allData.txt", sep="\t", as.is=T, header=T)[c(1,43)]
ped <- read.table("./pedigree/pedforQTLRel.txt", sep="\t", as.is=T, header=T)
names(geno.samples) <- "id"
geno.samples <- geno.samples + 0.1

load("./fmiss.RData")
fmiss <- fmiss[c(1,21)]
head(fmiss)


data.tmp <- merge(geno.samples, ail.ped, all.x=TRUE)
data.tmp$sire <- NULL; data.tmp$dam <- NULL
data.tmp$id <- data.tmp$id - 0.1
data.tmp <- merge(data.tmp, phenos, all.x=T)

data <- data.frame(ID_1=data.tmp$id, ID_2=data.tmp$id, missing=fmiss$meanMissing,
                   sex=data.tmp$sex,phenotype=data.tmp$act1.t)
headerline2 <- c(0,0,0,"D","P")
data <- rbind(headerline2, data)
head(data)

write.table(data, file="ailforPlink.sample", sep=" ", row.names=F, col.names=T, quote=F)


## plotting LD - one file method
library(ggplot2)
ld <- read.table("./dataFiles/chr19.ld",as.is=T, header=T)
ld$kb.dist <- with(ld, ld$BP_B-ld$BP_A)/1000  # gives distance in kb
#kb.distR <- signif(with(ld, ld$BP_B-ld$BP_A)/1000), digits=3)
meanr2 <- with(ld, tapply(ld$R^2, ld$kb.dist, FUN=mean, na.rm=T))
quantr95 <- with(ld, tapply(ld$R^2, ld$kb.dist, FUN=quantile, probs=0.95, na.rm=T))

ggplot(test)+
    geom_point(aes(x=distbins, y=meanr2))+
    geom_smooth(aes(x=distbins, y=quantr95), se=FALSE)+
    theme_bw()+
    xlab("Distance in 1kb bins")+
    ylab(expression(paste("Linkage disequilibrium (",r^2,")")))


## plotting LD script

setwd("/group/palmer-lab/AIL/GBS/dosage/")

files <- paste0("./chr", 1:19, ".ld")

ld.data <- function(file){
    temp <- read.table(file=file, as.is=T, header=T)[c(2,5,7)]
    temp$dist <- with(temp, temp$BP_B-temp$BP_A)/1000
    ld <- cbind.data.frame(temp$dist, temp$R)
    return(ld)
}

ld.info <- lapply(files, ld.data)
names(ld.info) <- paste0("chr", 1:19)
for(i in names(ld.info)){
    names(ld.info[[i]]) <- c("dist", "R")
}
linkageDf <- do.call(rbind.data.frame, ld.info)


meanr2 <- with(linkageDf, tapply(linkageDf$R^2, linkageDf$dist,
                                           FUN=mean, na.rm=T))
qt95 <- with(linkageDf, tapply(linkageDf$R^2, linkageDf$dist,
                                         FUN=quantile, probs=0.95, na.rm=T))
linkageD <- data.frame(meanr2, qt95)
linkageD$distbins <- as.numeric(rownames(linkageD))


library(ggplot2)
    pdf("LD_decay5.pdf", height=5, width=5.25)
    LDplot<- ggplot(linkageD)+
        geom_point(aes(x=distbins, y=meanr2))+
        geom_smooth(aes(x=distbins, y=qt95), se=FALSE)+
        theme_bw()+
        xlab("Distance (kB)")+
        ylab(expression(paste("Linkage disequilibrium (",r^2,")")))
    LDplot
    dev.off()


















