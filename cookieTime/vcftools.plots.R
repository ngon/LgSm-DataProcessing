chromosomes <- 1:19
cmds<-c()
for (chr in chromosomes){
    cmds[chr] <- paste0("vcftools --vcf ail.chr", chr, ".known.vcf --out chr", chr,
           " --site-depth, --site-mean-depth, --depth --geno-r2, --het")
}
write.table(cmds, "./getSeqStats.vcf.cmds", sep="\t", row.names=F, col.names=F, quote=F)

idepth <-c()
for (chr in chromosomes){
    idepth[[chr]] <- read.table(file=paste0("/group/palmer-lab/AIL/GBS/vars/vcf.stats/chr",chr,".idepth"), sep="\t", header=T, as.is=T)[2:3]
    idepth[[chr]] <- idepth[[chr]][do.call(order, idepth[[chr]]),]
}

library(ggplot2)

source("/group/palmer-lab/AIL/LgSm-DataProcessing/cookieTime/multiplot.R")
plotlist <- list()

for (chr in chromosomes) {

plotlist[[chr]] <- ggplot(data=idepth[[chr]], aes(x=idepth[[chr]]$N_SITES,
                                                y=idepth[[chr]]$MEAN_DEPTH))+

        geom_bar(stat="identity", fill="steelblue3") +
        xlab(paste0("Mean read depth per mouse on chr ", chr)) +
        ylab("Number of sites covered") +
        #scale_x_discrete(limits=(0:xmax), breaks=seq(0, xmax, 10)) +
        #scale_y_discrete(limits=(0:400), breaks=seq(0,400, 50))+
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=11),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))
}








pdf(file="idepth.chr.pdf", width=10, height=14)
multiplot(plotlist=plotlist, cols=2)
plot(idepth$N_SITES, idepth$MEAN_DEPTH, type='h',
     main=paste0("Read depth per mouse on chr ",chr),
     xlab="Number of sites", ylab="Mean read depth")
dev.off()


