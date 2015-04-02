
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
