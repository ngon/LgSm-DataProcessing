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

### IMPORT FILES ---------------------------------------------------------------
#
# hap <- read.table("./dataFiles/chr19.hap", header=F, as.is=T)
# names(hap) <- c("LG", "SM")
#
# legend <- read.table("./dataFiles/chr19.txt.legend.txt", header=F)[1]
# names(legend)<- c("snp", "ps", "ref", "alt")
#
# dos <- read.table("./dataFiles/chr19.filtered.dosage", header=F, as.is=T, nrows=200)[4:6]
# dosageFile <- "./dataFiles/chr19.filtered.dosage"
#
#
# ### FUNCTION: is.lgsm determines how many copies of the ref/alt alleles each mouse
# # (column) in the filtered.dosage file has at each locus.
#
# is.lgsm <- function(dosageFile){
#
#     genotypes <- read.table(dosageFile, header=F, as.is=T)[-c(1:3)]
#     genoClass <- c()
#
#     for (mouse in seq_along(genotypes)){
#         genoClass[[mouse]] <- cut(genotypes[[mouse]], breaks=c(0, 0.7, 1.3, 2),
#                          labels=c("R", "LS", "A"), dig.lab=4, right=TRUE,
#                          include.lowest=TRUE)
#     }
#     #genoClass <- do.call(what=cbind.data.frame, args=genoClass)
#     names(genoClass) <- 1:1830
#     genoClass <- lapply(genoClass, as.character)
#     return(genoClass)
# }
#
# testRun<- is.lgsm(dosageFile)


###############################################################################

is.lgsm <- function(chromosome){

    ######################## FIND WHICH STRAIN THE ALTERNATIVE ALLELE IS FROM
    legendFile <- paste0("/group/palmer-lab/AIL/knownSNPs/imputeHaplotypes/",
                      chromosome, ".txt")
    hapFile <- paste0("/group/palmer-lab/AIL/knownSNPs/imputeHaplotypes/",
                      chromosome, ".hap")

    legend <- read.table(legendFile, header=F, as.is=T)[1]
    haps <- read.table(hapFile, header=F, as.is=T)
    legend <- cbind(legend, haps)
    names(legend) <-  c("snp", "LG", "SM")
    rm(legendFile, hapFile, haps)

    which.alt <- c()
    which.ref <- c()
    for (i in legend[['LG']]){
        if (i == 1) {
            which.alt <- append(which.alt, "LL")
            which.ref <- append(which.ref, "SS")
        } else if (i == 0) {
            which.alt <- append(which.alt, "SS")
            which.ref <- append(which.ref, "LL")
        }
    }
    legend$which.alt <- as.character(which.alt)
    legend$which.ref <- as.character(which.ref)

    ################################################ CONVERT TO REF/NONREF
    dosageFile <- paste0("/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/",
                         chromosome, ".filtered.dosage")

    genotypes <- read.table(dosageFile, header=F, as.is=T)[-c(1:3)]

    #names(genotypes)[4:1833] <- paste0("V", 1:1830)
    genoClass <- c()

    for (mouse in seq_along(genotypes)){
        genoClass[[mouse]] <- cut(genotypes[[mouse]], breaks=c(0, 0.7, 1.3, 2),
                                  labels=c("R", "LS", "A"), dig.lab=4, right=TRUE,
                                  include.lowest=TRUE)
    }
    genoClass <- lapply(genoClass, as.character)
    #names(genoClass) <- 1:1830
    #################################
    snpnames <- read.table(dosageFile, header=F, as.is=T)[1]
    names(snpnames) <- c("snp")
    alt.allele <- legend[legend$snp %in% snpnames$snp,]
    ## figure out why these if statements only work in separate loops and not
    ## together under the same loop with 'else if (genoClass...)'
        for (mouse in seq_along(genoClass)) {
            for (i in seq_along(genoClass[[mouse]])){
                if (genoClass[[mouse]][i] == "A"){
                    genoClass[[mouse]][i] <- alt.allele$which.alt[i]
                }
            }
        }
    for (mouse in seq_along(genoClass)) {
        for (i in seq_along(genoClass[[mouse]])){
            if (genoClass[[mouse]][i] == "R"){
                genoClass[[mouse]][i] <- alt.allele$which.ref[i]
            }
        }
    }
    names(genoClass) <- 1:1830
    return(genoClass)
}








