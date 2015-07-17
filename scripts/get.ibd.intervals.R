### GET IBD INTERVALS ----------------------------------------------------------
### CODE TO MAKE LIST OF IBD INTERVALS BY CHROMOSOME FROM HEATHER LAWSON'S DATA


ibd.hl <- read.table("./dataFiles/LGSM.ibdRegions.hl.txt", sep="\t", header=T, as.is=T)
ibd.hl$Start <- round(ibd.hl$Start/1e6, digits=3)
ibd.hl$Stop <- round(ibd.hl$Stop/1e6, digits=4)

ibd <- ibd.hl[ibd.hl$State == "IBD",]
ibd <- split(ibd, ibd$Chr) ## note that chrs are not in numerical order
ibd[["chrX"]] <- NULL
rm(ibd.hl)
ibd.saved <- ibd
ibd <- ibd.saved

sumVariants <- c()
ibd.df <- c()
for (chr in seq_along(ibd)){
    for (i in 1:(nrow(ibd[[chr]])-1)){
        if(ibd[[chr]][["Stop"]][i] == ibd[[chr]][["Start"]][i+1]) {
            tmp <- ibd[[chr]][["Stop"]][i+1]
            ibd[[chr]][["Stop"]][i] = tmp
        }
        if(ibd[[chr]][["Stop"]][[i]] > ibd[[chr]][["Start"]][i+1]) {
            tmp <- ibd[[chr]][["Start"]][i]
            ibd[[chr]][["Start"]][[i+1]] <- tmp
        }
    }
    sumVariants[[chr]] <- tapply(ibd[[chr]][["Polymorphisms"]], ibd[[chr]][["Start"]], sum)
    ibd.df[[chr]] <- split(ibd[[chr]], ibd[[chr]][["Start"]])
}
names(ibd.df) <- names(ibd)


ibd.intervals <- vector("list", length=19)
for (chr in seq_along(ibd.df) ){
    for (j in seq_along(ibd.df[[chr]])) {
        ibd.intervals[[chr]][[j]] <- c(ibd.df[[chr]][[j]][nrow(ibd.df[[chr]][[j]]),2],
                                       ibd.df[[chr]][[j]][nrow(ibd.df[[chr]][[j]]),3])
    }
    ibd.intervals[[chr]] <- t(cbind.data.frame(ibd.intervals[[chr]]))
    row.names(ibd.intervals[[chr]]) <- 1:length(ibd.intervals[[chr]][,1])
    colnames(ibd.intervals[[chr]]) <- c("Start", "Stop")
    #intervals[[chr]] <- do.call(rbind.data.frame, intervals[[chr]])
}
ibd.intervals <- lapply(ibd.intervals, as.data.frame)
names(ibd.intervals) <- names(ibd)
rm(ibd, ibd.df, i, ibd.saved, j, chr)

# it's more convenient if the chromosomes are in numerical order
int <- list(chr1=ibd.intervals[["chr1"]], chr2=ibd.intervals[["chr2"]],
            chr3=ibd.intervals[["chr3"]], chr4=ibd.intervals[["chr4"]],
            chr5=ibd.intervals[["chr5"]], chr6=ibd.intervals[["chr6"]],
            chr7=ibd.intervals[["chr7"]], chr8=ibd.intervals[["chr8"]],
            chr9=ibd.intervals[["chr9"]], chr10=ibd.intervals[["chr10"]],
            chr11=ibd.intervals[["chr11"]], chr12=ibd.intervals[["chr12"]],
            chr13=ibd.intervals[["chr13"]], chr14=ibd.intervals[["chr14"]],
            chr15=ibd.intervals[["chr15"]], chr16=ibd.intervals[["chr16"]],
            chr17=ibd.intervals[["chr17"]], chr18=ibd.intervals[["chr18"]],
            chr19=ibd.intervals[["chr19"]])
ibd.intervals <- int
rm(int)
sumVariants <- list(sumVariants[[1]], sumVariants[[12]], sumVariants[[13]],
                    sumVariants[[14]], sumVariants[[15]], sumVariants[[16]],
                    sumVariants[[17]], sumVariants[[18]], sumVariants[[19]],
                    sumVariants[[2]], sumVariants[[3]], sumVariants[[4]],
                    sumVariants[[5]], sumVariants[[6]], sumVariants[[7]],
                    sumVariants[[8]], sumVariants[[9]], sumVariants[[10]], sumVariants[[11]])
names(sumVariants) <- paste0("chr", 1:19)

save(ibd.intervals, sumVariants, file="./dataFiles/ibd.intervals.RData")