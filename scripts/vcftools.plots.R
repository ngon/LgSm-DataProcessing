
### LOAD PLOTTING TOOLS ------------------------------------------------------
library(ggplot2)
source("/group/palmer-lab/AIL/LgSm-DataProcessing/scripts/multiplot.R")

### GENERATE VCFTOOLS SCRIPTS ------------------------------------------------
# chromosomes <- 1:19
# cmds<-c()
# for (chr in chromosomes){
#     cmds[chr] <- paste0("vcftools --vcf ail.chr", chr, ".known.vcf --out chr", chr,
#            " --site-depth, --site-mean-depth, --depth --geno-r2, --het")
# }
# write.table(cmds, "./getSeqStats.vcf.cmds", sep="\t", row.names=F, col.names=F, quote=F)

### PLOT COVERAGE PER MOUSE ---------------------------------------------------
idepth <-list()
for (chr in chromosomes){
    idepth[[chr]] <- read.table(file=paste0("/group/palmer-lab/AIL/GBS/vars/vcf.stats/chr",chr,".idepth"), sep="\t", header=T, as.is=T)[2:3]
    idepth[[chr]] <- idepth[[chr]][do.call(order, idepth[[chr]]),]
    mask <- apply(idepth[[chr]], 2, is.nan)
    is.na(idepth[[chr]][mask]) <- TRUE
  }

# make a list of ggplots for calling multiplot
pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/idepth.points.pdf",
    width=10, height=14)

plotlist <- list()
for (chr in seq_along(chromosomes)) {
plotlist[[chr]] <- ggplot(data=idepth[[chr]], aes(y=idepth[[chr]]$N_SITES,
                                                x=idepth[[chr]]$MEAN_DEPTH))+

        geom_point(stat="identity", position="identity", na.rm=T, color="orangered", alpha=0.5)+
        #geom_bar(stat="identity",position=position_dodge(width=0.8), fill="steelblue3") +
        xlab(paste0("Mean read depth per mouse on chr ", chr)) +
        ylab("Number of sites covered") +
        #scale_x_discrete(limits=(0:xmax), breaks=seq(0, xmax, 10)) +
        #scale_y_discrete(limits=(0:400), breaks=seq(0,400, 50))+
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10),
              axis.text.x=element_text(size=9),
              axis.text.y=element_text(size=9))
}


multiplot(plotlist=plotlist, cols=2)
dev.off()


### PLOT SUM AND MEAN COVERAGE PER SNP -----------------------------------------
snpdepth <-list()
for (chr in chromosomes){
    snpdepth[[chr]] <- read.table(file=paste0("/group/palmer-lab/AIL/GBS/vars/vcf.stats/chr", chr, ".ldepth"), sep="\t", header=T, as.is=T)[2:3]

    snpdepth[[chr]][3:4] <- read.table(file=paste0("/group/palmer-lab/AIL/GBS/vars/vcf.stats/chr",chr,".ldepth.mean"), sep="\t", header=T, as.is=T)[3:4]

    snpdepth[[chr]] <- snpdepth[[chr]][do.call(order, snpdepth[[chr]]),]
    mask <- apply(snpdepth[[chr]], 2, is.nan)
    is.na(snpdepth[[chr]][mask]) <- TRUE
    snpdepth[[chr]]$POS <- snpdepth[[chr]]$POS*1e-6
}
names(snpdepth)[1:19] <- chrnames[1:19]

### get locations of PstI cut sites -------------------------------------------

library(BSgenome.Mmusculus.UCSC.mm10)
ls("package:BSgenome.Mmusculus.UCSC.mm10")
genome <- BSgenome.Mmusculus.UCSC.mm10
psti = DNAString("CTGCAG")

chrnames <- paste0("chr", 1:19)
pstiSites <- list()
for (chr in chrnames){
    pstiSites[[chr]]<- matchPattern(pattern=psti, subject=DNAString(Mmusculus[[chr]]), max.mismatch=0, with.indels=F)
}
#rm(pstiSites,chrnames,genome,psti)
pstiCutPositions <- list()
for (chr in names(pstiSites)){
    pstiCutPositions[[chr]] <- start(pstiSites[[chr]])+4
}
pstiCutSites<- lapply(pstiCutPositions, data.frame, y=0)
for (chr in names(pstiCutSites)){
    names(pstiCutSites[[chr]])[1] <- "pos"
}
save(pstiCutSites, file="pstiCutSites_mm10.RData")

#
#     p<-    ggplot(data=snpdepth, aes(x=snpdepth$POS, y=log10(snpdepth$SUM_DEPTH)))+
#         geom_point(stat='identity', na.rm=T, color='blue', alpha=0.5)+
#         geom_point(data=snpdepth, aes(x=snpdepth$POS, y=log10(snpdepth$MEAN_DEPTH)),
#                        stat='identity', na.rm=T, color='orangered')

# make a list of ggplots for calling multiplot
plotlist <- list()
for (chr in seq_along(snpdepth)) {

m <- matrix(data=1:20, nrow=10, ncol=2)
pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/idepth.chr.pdf")
layout(m, widths=rep.int(1, ncol(m)), heights=rep.int(1, nrow(m)))

chr=18
plotlist[['chr18']] <-

   ggplot(data=snpdepth[[19]], aes(x=snpdepth[[19]]$POS,
                                  y=log10(snpdepth[[19]]$SUM_DEPTH))) +

        geom_point(stat="identity",
                   na.rm=T, color="blue", alpha=0.5)+

    geom_point(aes(x=snpdepth[[19]]$POS,
    y=log10(snpdepth[[19]]$MEAN_DEPTH)),
                   stat="identity",
                  na.rm=T, color="orangered", alpha=0.8) +

        xlab(paste0("mm10 Mb position on chr ", 19)) +
        ylab("Number of reads") +
        #scale_x_discrete(breaks=seq(0, max(snpdepth[[chr]]$POS)+2, 10)) +
        #scale_y_discrete(limits=(0:400), breaks=seq(0,400, 50))+
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=10),
              axis.text.x=element_text(size=9),
              axis.text.y=element_text(size=9))


}
# pdf(file="./test.pdf", width=10, height=14)
# multiplot(plotlist=plotlist, cols=2)
# dev.off()

pdf(file="/group/palmer-lab/AIL/LgSm-DataProcessing/figures/idepth.chr.pdf", width=10, height=14)
multiplot(plotlist=plotlist, cols=2)
dev.off()



