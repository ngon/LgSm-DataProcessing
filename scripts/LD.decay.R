## get LD decay plot
files <- paste0("./chr", 1:19, ".LD.10Mb.filtered")

ld.data <- function(file){
    temp <- read.table(file=file, as.is=T, header=T)[1:3] #[c(2,5,7)]
    temp$dist <- with(temp, temp$BP_B-temp$BP_A)/1000000
    ld <- cbind.data.frame(temp$dist, temp$R2)
    return(ld)
}

ld.info <- lapply(files, ld.data)
#names(ld.info) <- paste0("chr", 1:19)
for(i in seq_along(ld.info)){
    names(ld.info[[i]]) <- c("dist", "R2")
}
ld.info <- do.call(rbind.data.frame, ld.info)

meanr2 <- with(ld.info, tapply(ld.info$R2, ld.info$dist,FUN=mean, na.rm=T))
qt95 <- with(ld.info, tapply(ld.info$R2, ld.info$dist,FUN=quantile, probs=0.95, na.rm=T))
linkageD <- data.frame(meanr2, qt95)
rm(meanr2, qt95)
linkageD$distbins <- as.numeric(rownames(linkageD))


library(ggplot2)
LDplot<- ggplot(data=ld.info, aes(x=dist, y=R2))+
    geom_point(pch=".", color="ivory3")+
    geom_smooth(data=linkageD, aes(x=distbins, y=meanr2), se=FALSE, color="black")+
    geom_smooth(data=linkageD, aes(x=distbins, y=qt95), se=FALSE, color="slateblue")+
    theme_bw()+
    xlab("Distance (Mb)")+
    ylab(expression(paste("Linkage disequilibrium (",r^2,")")))

save(LDplot, file="./LDplotObject.RData")

pdf("/group/palmer-lab/AIL/LgSm-DataProcessing/figures/LD.10Mb_all.pdf",
    height=5, width=5.25)
LDplot
dev.off()


#
data <- read.table("./dataFiles/chr19.LD.10Mb", sep=" ", as.is=T, header=F)
names(data) <- c("BP_A", "BP_B", "R2")
data$dist <- with(data, data$BP_B-data$BP_A)/1000000
data <- cbind.data.frame(data$dist, data$R2)
names(data) <- c("dist", "R2")
ld.info <- data
meanr2 <- with(ld.info, tapply(ld.info$R2, ld.info$dist,FUN=mean, na.rm=T))
qt95 <- with(ld.info, tapply(ld.info$R2, ld.info$dist,FUN=quantile, probs=0.95, na.rm=T))
linkageD <- data.frame(meanr2, qt95)
rm(meanr2, qt95)
linkageD$distbins <- as.numeric(rownames(linkageD))

# library(ggplot2)
# LDplot<- ggplot(data=ld.info, aes(x=dist, y=R2))+
#     geom_point(pch=".", color="ivory3")+
#     theme_bw()+
#     xlab("Distance (Mb)")+
#     ylab(expression(paste("Linkage disequilibrium (",r^2,")")))
#
#     geom_smooth(data=linkageD, aes(x=distbins, y=meanr2), se=FALSE, color="black")+
#     geom_smooth(data=linkageD, aes(x=distbins, y=qt95), se=FALSE, color="slateblue")+
#
#
# save(LDplot, file="./LDplotObject.RData")
#
# pdf("/group/palmer-lab/AIL/LgSm-DataProcessing/figures/LD.10Mb_all.pdf",
#     height=5, width=5.25)
# LDplot
# dev.off()
