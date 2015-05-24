##### In progress: functions for summarizing genotype data ####

# IMPUTED VS EMPIRICAL SNPS ---------------------------------------------------

# # /AIL/GBS/preimpute/file.geno  <- empirical
# # read in a manageable part of chr 19 (smallest file)
# # cols: snpName, ref, alt, dosage...
# dos <- read.table("chr19.filtered.dosage", nrows=500)[1:500]
# names(dos)[1:3] <- c("snp", "ref", "alt")
#
# # /AIL/GBS/dosage/file.filtered.dosage <- imputed and empirical that passed QC
# # get names of first 500 empirical variants
# empirical <- read.table("../preimpute/ail.chr19.preimpute.geno", nrows=500)[2]
# names(empirical)[1] <- "snp"
#
# # get separate df for empirical and imputed snps
# imp <- dos[!dos$snp %in% emp$snp,]
# emp <- dos[dos$snp %in% emp$snp,]
#
# # write on CRI, push to github and load here for testing
# imp <- read.table("imputed_19.txt", sep="\t", header=T)
# emp <- read.table("empirical_19.txt", sep="\t", header=T)

# GET NAMES OF EMPIRICAL AND IMPUTED SNPS -----------------------------------
# Run R from /group/palmer-lab/AIL/GBS/genoSummaries

#get SNP names...the file is huge and it takes a long time ( >60 min)
    allSnps <- read.table(file="../dosage/chrAll.filtered.dosage",
                           sep=" ",as.is=TRUE)[1]
    names(allSnps)[1] <- "snps"
    write.table(allSnps, file="./quality.snp.list.txt", sep="/t", col.names=TRUE,
                row.names=FALSE, quote=FALSE) # 3479041 SNPs

    # get preimpute filenames except for chrX (idx 20), which is empty
    preimpFiles <- list.files(path="../preimpute/", pattern="*.preimpute.geno")
    preimpFiles <- preimpFiles[-20]

    # extract SNP names from each file. this takes quite a long time. (20 min)
    get.emp.positions <- function(file) {read.table(file, sep="\t", header=F)[3]}
    empSnps <- lapply(file.path("../preimpute/", preimpFiles), get.emp.positions)
    empSnps <- do.call(what=rbind.data.frame, args=empSnps) # 316,013 SNPs
    names(empSnps)[1] <- "snps"
    write.table(empSnps, file="./empirical.snp.list.txt", sep="/t",quote=FALSE,
            row.names=FALSE, col.names=TRUE)

    # get separate lists of empirical and imputed SNPs
    emp <- which(allSnps$snps %in% empSnps$snps) # 312516
    imp <- which(!allSnps$snps %in% empSnps$snps) # 3166525

empRows=list()
for(file in filenames){
    empRows[[file]] <- which((read.table(file, header=F, sep="\t",as.is=T)[2])$V2 %in% emp$ps)
}
#  MAF and HET ---------------------------------------------------------------
# is there anything <16% (copying error rate in imputation)

# list the necessary file names
# chromosomes <- paste0("chr", 1:19)
# filenames <- list()
# for (i in chromosomes){
#     filenames[i] <- paste0("../dosage/", i, ".filtered.dosage")
# }

# load positions of empirical snps (a list called empRows)
chromosomes <- c("chr8", "chr19")
filenames <- c("./chr8.txt", "./chr19.txt")
names(filenames) <- c("chr8", "chr19")

load("empSnpRows.Rdata")

# GET MAF AND HET ----------------------------------------------------------------------

### get.empirical.data ----------------------------------------------------
get.empirical.data <- function(file, e=emp.rows, x=maf, y=het.snp){
    print("Getting data on empirical SNPs...")

    #chrname <- substr(basename(file), 1, nchar(basename(file))-16)
    chrname <- substr(basename(file), 1, nchar(basename(file))-4)
    empR <- e[[chrname]]
    #print(chrname)
    #print (length(empR))
    #emp.maf <- data.frame(table(sort(x[seq_along(x) %in% empR])))
    emp.maf <- sort(x[seq_along(x) %in% empR])
    #names(emp.maf)[1] <- "emp.maf"
    empmaf.hist <- hist(emp.maf[1], breaks=seq(0,0.5, 0.05), plot=F,
                        include.lowest=T)
    #emp.het <- data.frame(table(sort(y[seq_along(y) %in% empR])))
    emp.het <- sort(y[seq_along(y) %in% empR])
    #names(emp.het)[1] <- "emp.het.snp"
    emphet.hist <- hist(emp.het[1], breaks=seq(0,0.5, 0.05), plot=F,
                        include.lowest=T)


   # list(emp.maf=emp.maf, emp.het=emp.het)
    list(empmaf.hist=empmaf.hist, emphet.hist=emphet.hist)
}

### maf and het --------------------------------------------------------------
maf.and.het <- function(file, upper=1.2, lower=0.8, emp.rows=NULL) {
#     chrname <- names(file)
    print("Getting genotype data...")
    geno <- read.table(file, header=F, as.is=T)[-c(1:4)] # change 4 to a 3

    ########## MAF
    print("Calculating minor allele frequencies...")
    freq <- cbind( (rowMeans(geno, na.rm = T)/2),
                1-(rowMeans(geno, na.rm=T)/2)  )
    maf <- round(apply(freq, 1, min), digits=2) ####
    maf.hist <- hist(maf, breaks=seq(0,0.5, 0.05), plot=F,
                     include.lowest=T) #####

    ########## HET
    print("Calculating average heterozygosity...")
    # get heterozygosity by sample and by SNP
    hetsites <- geno <= upper & geno >= lower
    het.mouse <- apply(hetsites, 2, mean)
    het.snp <- round(apply(hetsites, 1, mean), digits=2) ####
    hetsnp.hist <- hist(het.snp, breaks=seq(0,0.5, 0.05), plot=F,
                        include.lowest=T) #####

    if(!is.null(emp.rows)){
        empStats <- get.empirical.data(file=file, e=emp.rows, x=maf, y=het.snp)

        #maf <- data.frame(table(sort(maf)))
        #names(maf)[1] <- "maf"
        #het.snp <- data.frame(table(sort(het.snp)))
        #names(het.snp)[1] <- "het.snp"

        result<- list(maf.hist, hetsnp.hist, het.mouse,
                  empStats["empmaf.hist"], empStats["emphet.hist"])

    } else {
        #maf <- data.frame(table(sort(maf)))
        #names(maf)[1] <- "maf"
        #het.snp <- data.frame(table(sort(het.snp)))
        #names(het.snp)[1] <- "het.snp"

        result<-  list(maf, het.snp, het.mouse)
    }
    return(result)

}

## plot maf  ------------------------------------------------------
data <- lapply(filenames, maf.and.het, emp.rows=empRows)

# files <- c("./chr8.txt", "./chr19.txt")
# names(files)[1:2]<- c("chr8", "chr19")
# top level test = chr
# bottom level = maf, het.snp, het.mouse (vector), emp.maf, emp.het


pdf(file="./mafc.pdf", height=5, width=10, colormodel="cmyk")
par(mfrow=c(1,2))
for(chr in chromosomes){
    #all <- data.frame(data[[chr]][1])
    all <- data[[chr]][1]
    #names(all)[1:2] <-c("maf", "Freq")
    #emp <- data.frame(data[[chr]][4])
    emp <- data[[chr]][4]
    #names(emp)[1:2] <-c("maf", "Freq")
    plot.het.maf(all, emp)

}
dev.off()

plot.het.maf <- function(all, emp) {
   plot(x=all[[1]], xlab=paste0("Minor allele frequency on ", chr),
             ylab="Number of SNPs", main=NULL, cex.lab=0.9, cex.axis=0.8,
             col=rgb(255,193,37, maxColorValue = 255))

   plot(x=emp[[1]]$empmaf.hist, col=rgb(255,69,0,128, maxColorValue = 255), add=T)

   legend(x=0, y=max(all[[1]]$counts), legend=c("All", "GBS"), bty="n",
               cex=0.8, fill=c(rgb(255,193,37, maxColorValue=255),
                               rgb(255,69,0,150, maxColorValue=255)))
    box()

    }


    mafplots[[chr]] <- ggplot(data=all, aes(x=maf, y=Freq/sum(Freq))) +
            geom_bar(stat="identity", fill="goldenrod1", color="black") +
            geom_bar(data=emp, aes(x=maf, y=Freq/sum(Freq)),
                     stat="identity", fill="orangered", alpha=0.5, color="black") +
            xlab(paste0("Minor allele frequency on ", chr)) +
            ylab("Proportion of SNPs") +
            theme_bw() +
            theme(axis.title.x=element_text(size=10),
                  axis.title.y=element_text(size=11),
                  axis.text.x=element_text(size=8),
                  axis.text.y=element_text(size=10))

    print("Plotting SNP heterozygosity for each chromosome....")
    all <- data.frame(data[[chr]][2])
    names(all)[1:2] <-c("het", "Freq")
    emp <- data.frame(data[[chr]][5])
    names(emp)[1:2] <-c("het", "Freq")

    hetplots[[chr]] <- ggplot(data=all, aes(x=het, y=Freq/sum(Freq))) +
        geom_bar(stat="identity", fill="goldenrod1", color="black") +
        geom_bar(data=emp, aes(x=het, y=Freq/sum(Freq)),
                 stat="identity", fill="orangered", alpha=0.5, color="black") +
        xlab(paste0("Average heterozygosity on ", chr)) +
        ylab("Proportion of SNPs") +
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=11),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))

    print("Plotting mouse heterozygosity for each chromosome....")
    all <- as.data.frame(unlist(data[[chr]][3]))
    names(all)[1] <- "het"

    mouseplots[[chr]] <- ggplot(data=all, aes(x=het)) +
        geom_histogram(fill="goldenrod1", color="black", binwidth=0.1) +
        xlab(paste0("Average heterozygosity on ", chr)) +
        ylab("Proportion of mice") +
        theme_bw() +
        theme(axis.title.x=element_text(size=10),
              axis.title.y=element_text(size=11),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))

}



pdf(file="./mafChromosomes.pdf", height=12, width=10, colormodel="cmyk")
multiplot(plotlist=mafplots, cols=4)
dev.off()

pdf(file="./hetChromosomes.pdf", height=12, width=10, colormodel="cmyk")
multiplot(plotlist=hetplots, cols=4)
dev.off()

pdf(file="./hetmouseChromosomes.pdf", height=10, width=8, colormodel="cmyk")
multiplot(plotlist=mouseplots, cols=5)
dev.off()


# GRM RELATEDNESS ------------------------------------------------------
# first centralize rows of the pxn matrix where p=snps, n=mice
# a. find mean dosage for each SNP
# b. subtract the mean dosage from every mouse's dose
# c. GRM = ( t(Xc) %*% Xc ) / p, where Xc is the centralized version of the data
#    file and P is the number of SNPs. this is the var/covar matrix.

geno <- read.table("chr19.txt")[-c(1:4)]
meanDosage <- rowMeans(geno, na.rm = T)

# the matrix has to be transposed inside the scale() function because scale()
# does column-wise centering.
centerGeno <- (scale(t(as.matrix(geno)), center=meanDosage, scale=FALSE))

# when you simply use cov or var(matEmp), you get slightly different results
GRM <- centerGeno%*%(t(centerGeno)) / nrow(centerGeno)


# LD DECAY -------------------------------------------------------------------
# transpose GRM correl matrix to get cor(r) between snps
# pairwise LD: (cor(Xrow1, Xrow2...))^2
# AIL/knownSNPs/imputeHaplotypes .hap files, first snp is LG and 2nd is SM.
# DOSAGE = always the number of alternative alleles, so you need to match back to
# LG and SM using the genos in the haplo files.
# .legend files give row names (snps)

# compute LD
LD <- (t(centerGeno)%*%centerGeno/ncol(centerGeno))**2

# get snp names
snps <- read.table("chr19.txt")[2]
snps <- as.list(substring(snps$V2, first=11))
rownames(LD)[1:500] <- snps
colnames(LD)[1:500] <- snps


# to plot LD decay, i want a df with SNP pairs in the first two cols, positions,
# and the correlation between them
# distance = site.i - site.j

flattenMatrix <- function(m) {
    data.frame(snp1=as.double(as.character(rep(row.names(m),ncol(m)))),
                snp2=as.double(as.character(rep(colnames(m),each=nrow(m)))),
                rsq=as.vector(m))
}

test <- flattenMatrix(LD)
test$dist <- abs(test$snp2 - test$snp1)

row_bp <- unique(test$snp1)
col_bp <- unique(test$snp2)



