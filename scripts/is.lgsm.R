### June 15, 2015
# PURPOSE: Determine whether alternative allele in .filtered.dosage files
# is LG or SM by linking to info in haplotype files.
## THIS SCRIPT HAS BEEN TESTED ON CRI. EXECUTE ONE CHROMOSOME AT A TIME USING:
## > is.lgsm(chromosome="chr19")
## OR ALL AT ONCE WITH: > haplos <- lapply(chrNameVector, is.lgsm)

### DESCRIPTION OF DATA --------------------------------------------------------
# /AIL/knownSNPs/imputeHaplotypes/ contains .hap and .legend files.
# .hap files have 2 cols - col 1 is LG and col 2 is SM. values in these cols
# are either 0 for the reference allele, or 1 for the alternative allele.
# Legend files list the reference and alternative alleles' nucleotide bases.
# Dosage in .filtered.dosage files corresponds to the number of alternative
# alleles each individual has at the locus.

# filtered dosage files are rows of SNPs. first 3 cols are same as .legend file
# write function to act on each element in a vector
# apply function to each vector in a data frame


# ### FUNCTION: is.lgsm determines how many copies of the ref/alt alleles each mouse
# (column) in the filtered.dosage file has at each locus. it matches ref and alt to
# LG and SM, then creates a coded list of LG (2), SM (0), or HET (1) alleles for each
# genotyped sample to be used for plotting and analysis.

is.lgsm <- function(chromosome,
                    legendPathPrefix = "/group/palmer-lab/AIL/knownSNPs/imputeHaplotypes/",
                    haploPathPrefix = "/group/palmer-lab/AIL/knownSNPs/imputeHaplotypes/",
                    filtDosPathPrefix = "/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/")  {


    legendFile      <- paste0(legendPathPrefix, chromosome, ".txt")
    hapFile         <- paste0(haploPathPrefix, chromosome, ".hap")

    # FIND WHICH STRAIN THE ALTERNATIVE ALLELE IS FROM
    legend          <- read.table(legendFile, header=F, as.is=T)[1]
    haps            <- read.table(hapFile, header=F, as.is=T)
    legend          <- cbind(legend, haps)
    names(legend)   <-  c("snp", "LG", "SM")
    rm(legendFile, hapFile, haps)

    which.alt <- c()
    which.ref <- c()
    for (i in legend[['LG']]){
        if (i == 1) {
            which.alt <- append(which.alt, 2) # 2 = 2 copies of LG
            which.ref <- append(which.ref, 0) # 0 = 0 copies of LG (2 of SM)
        } else if (i == 0) {
            which.alt <- append(which.alt, 0)
            which.ref <- append(which.ref, 2)
        }
    }
    legend$which.alt <- as.character(which.alt)
    legend$which.ref <- as.character(which.ref)

    # DETERMINE WHETHER EACH ALLELE IS REF OR ALT
    dosageFile <- paste0(filtDosPathPrefix, chromosome, ".filtered.dosage")

    genotypes <- read.table(dosageFile, header=F, as.is=T)[-c(1:3)]

    genoClass <- c()

    for (mouse in seq_along(genotypes)){
        genoClass[[mouse]] <- cut(genotypes[[mouse]], breaks=c(0, 0.8, 1.2, 2),
                                  labels=c("R", 1, "A"), dig.lab=4, right=TRUE,
                                  include.lowest=TRUE) # 1 = LgSm het
    }
    genoClass <- lapply(genoClass, as.character)

    # REPLACE 'A' AND 'R' WITH LG/SM CODED VALUES
    snpnames <- read.table(dosageFile, header=F, as.is=T)[1]
    names(snpnames) <- c("snp")
    alt.allele <- legend[legend$snp %in% snpnames$snp,]
    ## figure out why these if statements only work in separate loops and not
    ## together under the same loop with 'else if (genoClass...)'
        for (mouse in seq_along(genoClass)) {
            for (i in seq_along(genoClass[[mouse]])){
                if (genoClass[[mouse]][i] == "A"){
                    genoClass[[mouse]][i] <- alt.allele$which.alt[i]
                } # if
            } # for i
        } # for mouse
    for (mouse in seq_along(genoClass)) {
        for (i in seq_along(genoClass[[mouse]])){
            if (genoClass[[mouse]][i] == "R"){
                genoClass[[mouse]][i] <- alt.allele$which.ref[i]
            } # if
        } # for i
    } # for mouse
    names(genoClass) <- 1:1830
    genoClass <- lapply(genoClass, as.integer)
    return(genoClass)
}

############ RUNNING isLgorSm.R ################################################
# ids <- read.table("/group/palmer-lab/AIL/LgSm-DataProcessing/genotyped.samples.txt", as.is=T)$V1

# the commands have to be elaborated individually because i want each chr in its
# own separate file. a combined file would be too large for R to keep in memory.

# chr1.hap <- is.lgsm(chromosome="chr1")
# names(chr1.hap) <- ids
# save(chr1.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr1.hap.RData")
# rm(chr1.hap)
#
# chr2.hap <- is.lgsm(chromosome="chr2")
# names(chr2.hap) <- ids
# save(chr2.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr2.hap.RData")
# rm(chr2.hap)
#
# chr3.hap <- is.lgsm(chromosome="chr3")
# names(chr3.hap) <- ids
# save(chr3.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr3.hap.RData")
# rm(chr3.hap)
#
# chr4.hap <- is.lgsm(chromosome="chr4")
# names(chr4.hap) <- ids
# save(chr4.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr4.hap.RData")
# rm(chr4.hap)
#
# chr5.hap <- is.lgsm(chromosome="chr5")
# names(chr5.hap) <- ids
# save(chr5.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr5.hap.RData")
# rm(chr5.hap)
#
# chr6.hap <- is.lgsm(chromosome="chr6")
# names(chr6.hap) <- ids
# save(chr6.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr6.hap.RData")
# rm(chr6.hap)
#
# chr7.hap <- is.lgsm(chromosome="chr7")
# names(chr7.hap) <- ids
# save(chr7.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr7.hap.RData")
# rm(chr7.hap)
#
# chr8.hap <- is.lgsm(chromosome="chr8")
# names(chr8.hap) <- ids
# save(chr8.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr8.hap.RData")
# rm(chr8.hap)
#
# chr9.hap <- is.lgsm(chromosome="chr9")
# names(chr9.hap) <- ids
# save(chr9.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr9.hap.RData")
# rm(chr9.hap)
#
# chr10.hap <- is.lgsm(chromosome="chr10")
# names(chr10.hap) <- ids
# save(chr10.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr10.hap.RData")
# rm(chr10.hap)
#
# chr11.hap <- is.lgsm(chromosome="chr11")
# names(chr11.hap) <- ids
# save(chr11.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr11.hap.RData")
# rm(chr11.hap)
#
# chr12.hap <- is.lgsm(chromosome="chr12")
# names(chr12.hap) <- ids
# save(chr12.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr12.hap.RData")
# rm(chr12.hap)
#
# chr13.hap <- is.lgsm(chromosome="chr13")
# names(chr13.hap) <- ids
# save(chr13.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr13.hap.RData")
# rm(chr13.hap)
#
# chr14.hap <- is.lgsm(chromosome="chr14")
# names(chr14.hap) <- ids
# save(chr14.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr14.hap.RData")
# rm(chr14.hap)
#
# chr15.hap <- is.lgsm(chromosome="chr15")
# names(chr15.hap) <- ids
# save(chr15.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr15.hap.RData")
# rm(chr15.hap)
#
# chr16.hap <- is.lgsm(chromosome="chr16")
# names(chr16.hap) <- ids
# save(chr16.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr16.hap.RData")
# rm(chr16.hap)
#
# chr17.hap <- is.lgsm(chromosome="chr17")
# names(chr17.hap) <- ids
# save(chr17.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr17.hap.RData")
# rm(chr17.hap)
#
# chr18.hap <- is.lgsm(chromosome="chr18")
# names(chr18.hap) <- ids
# save(chr18.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr18.hap.RData")
# rm(chr18.hap)
#
# chr19.hap <- is.lgsm(chromosome="chr19")
# names(chr19.hap) <- ids
# save(chr19.hap, file="/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/chr19.hap.RData")
# rm(chr19.hap)
