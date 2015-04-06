#### plot.SNP.density ######
#### a function for plotting SNP density for each chromosome
#### Run this in /group/palmer-lab/AIL/GBS/genoSummaries


# get preimpute filenames except for chrX (idx 20), which is empty
preimpFiles <- list.files(path="../preimpute/", pattern="*.preimpute.geno")
preimpFiles <- preimpFiles[-20]
# extract row names from each file. the resulting list is 1.2 Mb
get.emp.rownames <- function(file) {read.table(file, sep="\t", header=F)[3]}
empSnps <- lapply(file.path("../preimpute/", preimpFiles), get.emp.rownames)
# sort list
goodOrder <- paste0("chr", 1:19)
empSnp <- list()
for(chr in goodOrder) {
    empSnp[[chr]] <- empSnps[[chr]]
}
emp <- do.call(rbind.data.frame, empSnps)
names(emp)[1] <- "ps"
save(empSnps, emp, file="empSnpPositions.RData")

chrlens <- trunc((read.table("chrLengths_mm10.txt")[2])/1e6)
save(emp, chrlens, file="empSnps.and.chrlens.RData")

##############################################################################
# to do: read in each chromosome at once, get the information you need and then
# discard from memory.
# all you need is ps col from each chrInfo file. immediately divide all the ps
# by 1 Mb (or whatever bin size you want) and store the results AS INTEGERS.
# then: chr1bins<- tapply(sort(unique(chrInfoTable)), chrInfoTable, count)
# this will give you the number of SNPs in each 1 Mb bin. save this into a list
# use hist(chr1bins, plot=FALSE)
# then match by row name to list of empirical snps

######## SETUP -----------------------------------------------------------------

library("ggplot2")
# load list of empirical SNP positions (emp) and chromosome lengths (chrlens)
load("empSnps.and.chrlens.RData")

snp.density <- function(filename, empRows=FALSE) {
    chrPos <- trunc((read.table(filename, sep="\t", header=F)[2])/1e6)
    allDensity <- data.frame(table(sort(chrPos$V2)))
    names(allDensity)[1] <- "all.Mb"

    if(empRows==TRUE){
        is.emp <- which((read.table(filename, sep="\t", header=F, as.is=T)[2])$V2 %in% emp$ps)
        empDensity <- data.frame(table(sort(chrPos[row.names(chrPos) %in% is.emp,])))
        names(empDensity)[1] <- "emp.Mb"
        list(allDensity, empDensity)
    } else {
        list(allDensity)
    }

    }

####### COMMANDS  --------------------------------------------------------------
chromosomes <- paste0("chr", 1:19)
#chromosomes <- paste0("chr", 17:19)
filenames <- list()
#for (i in chromosomes) { filenames[i] <-paste0(i, "_info_short")}
     for (i in chromosomes){
         filenames[i] <- paste0("../dosage/", i, ".filtered.snpinfo")
     }

snpsPerMb <- list()
for(file in filenames) {
    snpsPerMb[[file]] <- snp.density(file, empRows=TRUE)
    }
names(snpsPerMb)[1:19] <- paste0("chr", 1:19)
#names(snpsPerMb)[1:3] <- paste0("chr", 17:19)

####### PLOTS ------------------------------------------------------------------
chrlens.test <- chrlens[17:19,]
chrdata <- snpsPerMb[["chr19"]]

plotlist <- list()

for (chr in seq_along((snpsPerMb))) {
    chrdata <- snpsPerMb[[chr]]
    names(chrdata)[1:2] <- c("all", "emp")
    #xmax <- chrlens[chr,]
    xmax <- chrlens.test[chr]

plotlist[[chr]] <- ggplot(data=chrdata[[1]], aes(x=all.Mb, y=Freq)) +
        geom_bar(stat="identity", fill="goldenrod2") +
        geom_bar(data=chrdata[[2]], aes(x=emp.Mb, y=Freq),
                       stat="identity", fill="steelblue3") +
        #ggtitle("Genome-wide SNP density") +
        xlab(paste0("Position (Mb) on chromosome ", chr)) +
        ylab("SNP density") +
        scale_x_discrete(limits=(0:xmax), breaks=seq(0, xmax, 10)) +
        #scale_y_discrete(limits=(0:400), breaks=seq(0,400, 50))+
        theme_bw() +
        theme(axis.title.x=element_text(size=11),
          axis.title.y=element_text(size=11),
          axis.text.x=element_text(size=10),
          axis.text.y=element_text(size=10))
}

pdf(file="./snpDensityPlots.pdf", height=14, width=10)
multiplot(plotlist=plotlist, cols=2)
dev.off()







