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

get SNP names...the file is huge and it takes a long time ( >60 min)
    allSnps <- read.table(file="../dosage/chrAll.filtered.dosage",
                           sep=" ",as.is=TRUE)[1]
    names(allSnps)[1] <- "snps"
    write.table(allSnps, file="./quality.snp.list.txt", sep="/t", col.names=TRUE,
                row.names=FALSE, quote=FALSE) # 3479041 SNPs

    # get preimpute filenames except for chrX (idx 20), which is empty
    preimpFiles <- list.files(path="../preimpute/", pattern="*.preimpute.geno")
    preimpFiles <- preimpFiles[-20]

    # extract SNP names from each file. this takes quite a long time. (20 min)
    get.emp.rownames <- function(file) {read.table(file, sep="\t", header=F)[3]}
    empSnps <- lapply(file.path("../preimpute/", preimpFiles), get.emp.rownames)
    empSnps <- do.call(what=rbind.data.frame, args=empSnps) # 316,013 SNPs
    names(empSnps)[1] <- "snps"
    write.table(empSnps, file="./empirical.snp.list.txt", sep="/t",quote=FALSE,
            row.names=FALSE, col.names=TRUE)

    # get separate lists of empirical and imputed SNPs
    emp <- which(allSnps$snps %in% empSnps$snps) # 312516
    imp <- which(!allSnps$snps %in% empSnps$snps) # 3166525


#  MAF and HET ---------------------------------------------------------------
# is there anything <16% (copying error rate in imputation)

# geno is the filtered.dosage file
maf.and.het <- function(genofile, empList=NULL, color=rgb(0.1,0.1,0.1,0.5),
                    color2=rgb(1,0,0,0.5)) {

    print("Getting genotype data...")
    geno <- read.table(genofile, header=F, as.is=T, nrows =3)[-c(1:3)]
    classes <- lapply(geno, class)
    geno <- as.matrix(read.table(genofile, header=F, colClasses=classes))

    if(!is.null(empList)){
        emp <- read.table(file=empList, sep="\t", header=T)[1]
    }

    print("Calculating minor allele frequencies...")
    freq <- cbind((rowMeans(geno, na.rm = T)/2),
                  1-(rowMeans(geno, na.rm=T)/2))
    maf <- apply(freq, 1, min)
    return(maf)

    print("Plotting MAF histogram....")
    png(file="MAFplot.png", width=425, height=375, units="px")

    if(is.null(empList)) { # Standard MAF histogram
         hist(maf, breaks=seq(0, 0.5, by=0.01), col=color,
                    main = "Minor allele frequency distribution",
                    xlab = "Minor allele frequency",
                    ylab = "Count")
         box()
    }

    else { # Draw MAF histogram stratified by empirical vs. imputed SNPs
         hist(maf[which(seq_along(maf) %in% emp[1])],
                    breaks=seq(0, 0.5, by=0.01), col=color,
                    main = "Minor allele frequency distribution",
                    sub  = "Empirical vs. imputed SNPs",
                    xlab = "Minor allele frequency",
                    ylab = "Count")
        hist(maf[which(!seq_along(maf) %in% emp[1])],col=color2, add=TRUE)
        box()
           }
    dev.off()

    print("Calculating average heterozygosity...")
    # get heterozygosity by sample and by SNP
    hetsites <- geno <= upper & geno >= lower
    het.mouse <- apply(hetsites, 2, mean)
    het.snp <- apply(hetsites, 1, mean)
    list(het.mouse, het.snp)

    print("Making heterozygosity plots...")
    png(file="het.mouse.png", width=425, height=375, units="px")

    # het.mouse histogram
    hist(het.mouse, breaks=30, col=color,
             main = "Average heterozygosity per mouse",
             xlab = "Average heterozygosity",
             ylab = "Number of mice")
    dev.off()

    png(file="het.snp.png", width=425, height=375, units="px")

    if(is.null(empList)) { # Standard het.snp histogram
        hist(het.snp, breaks=30, col=color,
             main = "Average heterozygosity per SNP",
             xlab = "Average heterozygosity",
             ylab = "Number of SNPs")
        box()
    }

    else { # Het.snp histogram stratified by empirical vs. imputed SNPs
        hist(het.snp[which(seq_along(het.snp) %in% emp[1])],
             breaks= 30, col=color,
             main = "Average heterozygosity per SNP",
             sub  = "Empirical vs. imputed SNPs",
             xlab = "Average heterozygosity",
             ylab = "Number of SNPs")
        hist(maf[which(!seq_along(maf) %in% emp[1])],col=color2, add=TRUE)
        box()
    }

    dev.off()
    }

# HETEROZYGOSITY BY SAMPLE ----------------------------------------------------
# for each row (snp) calculate the proportion of cols with dosage 0.8-1.2.

heterozygosity <- function(geno, lower=0.8, upper=1.2){

    df <- as.matrix(geno[-c(1:3)])
    hetsites <- df <= upper & df >= lower

    het.mouse <- apply(hetsites, 2, mean)
    het.snp <- apply(hetsites, 1, mean)

    list(het.mouse, het.snp)
    }

# plot HET
hist(ehet.mouse, breaks=30, main="Average heterozygosity per mouse",
     xlab="Average heterozygosity")

# HETEROZYGOSITY BY SNP -------------------------------------------------------
# for each row (snp) calculate the proportion of cols with dosage 0.8-1.2.

het.snp <- function(geno, lower=0.8, upper=1.2){
    df <- geno[-c(1:3)]
    hets <- vector(length=nrow(df))
    for(i in 1:nrow(df)) {
        hets[i] <- (sum(df[i,] >= lower & df[i,] <= upper))/(ncol(df))
    }
    return(hets)
}
# plot HET
hist(ehet.snp, breaks=30, main="Average heterozygosity per SNP",
     xlab="Average heterozygosity")

# LOCAL NUCLEOTIDE DIVERSITY ---------------------------------------------------
# how many SNPs in each n Kb window - write a fxn where the user can specify n
# code SNPs by strain origin
# might also be interesting to make a 3D plot including association pvals



# PED VS GRM RELATEDNESS ------------------------------------------------------
# 1. Get GRM
# first centralize rows of the pxn matrix where p=snps, n=mice
### eg subtract the mean dosage for SNPi from every mouse's dose
# GRM = ( t(Xc) %*% Xc ) / p, where Xc is the centralized version of the data file and P
# is the number of SNPs. this is the var/covar matrix.
# to get pairwise LD, you want the (cor(Xrow1, Xrow2...))^2

meanDosage <- rowMeans(emp[-c(1:3)], na.rm = T)

# the matrix has to be transposed inside the scale() function because scale()
# does column-wise centering.
matEmp <- t(scale(t(as.matrix(emp[-c(1:3)])), center=meanDosage, scale=FALSE))

empGRM <- (t(matEmp)%*%(matEmp)) / nrow(matEmp)
# when you simply use cov or var(matEmp), you get slightly different results


# 2. Compare to ped


# LD DECAY -------------------------------------------------------------------
# transpose GRM correl matrix to get cor(r) between snps
# AIL/knownSNPs/imputeHaplotypes
# .hap files, first snp is LG and 2nd is SM.
# DOSAGE = always the number of alternative alleles, so you need to match back to
# LG and SM using the genos in the haplo files.
# .legend files give row names (snps)

empLD <- t(cor(t(empGRM)))^2


# to plot LD decay, i want a df with SNP pairs in the first two cols, positions,
# and the correlation between them
# distance = site.i - site.j
# order df by distance
# row_bp <- unique(df$site.i); col_bp <- unique(df$site.j)







