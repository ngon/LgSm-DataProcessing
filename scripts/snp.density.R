#### Purpose: get SNP density data and plot it for each chromosome

##### GET LIST OF EMPIRICAL SNPS
## This is how the list of empirical SNPs was obtained. The "full" set of
## SNPs is from the first pass genotypes called using Cheverud LG/SM known SNPs.
# get preimpute filenames except for chrX (idx 20), which is empty
# preimpFiles <- list.files(path="../preimpute/", pattern="*.preimpute.geno")
# preimpFiles <- preimpFiles[-20]
# # extract row names from each file. the resulting list is 1.2 Mb
# get.emp.rownames <- function(file) {read.table(file, sep="\t", header=F)[3]}
# empSnps <- lapply(file.path("../preimpute/", preimpFiles), get.emp.rownames)
# # sort list
# goodOrder <- paste0("chr", 1:19)
# empSnp <- list()
# for(chr in goodOrder) {
#     empSnp[[chr]] <- empSnps[[chr]]
# }
# emp <- do.call(rbind.data.frame, empSnps)
# names(emp)[1] <- "ps"
# save(empSnps, emp, file="empSnpPositions.RData")

###################################################
#### FUNCTION: snp.density
#### Requires a list of lists for each chromosome with positions of empSNPs.

snp.density <- function(chr, empRows=FALSE, empPositions=NULL) {

    filename <- paste0("/group/palmer-lab/AIL/GBS/dosage/", chr, ".filtered.snpinfo")
    chrPos <- trunc((read.table(filename, sep="\t", header=F)[2])/1e6)
    allDensity <- data.frame(table(sort(chrPos$V2)))
    names(allDensity)[1] <- "all.Mb"

    if(empRows==TRUE){
        is.emp <- which((read.table(filename, sep="\t", header=F)[2])$V2 %in% empPositions[[1]])
        empDensity <- data.frame(table(sort(chrPos[row.names(chrPos) %in% is.emp,])))
        names(empDensity)[1] <- "emp.Mb"
        list(allDensity, empDensity)
    } else {
        list(allDensity)
    }

}
### commands to run snp.density
# load("./empSnpsPosByChr.RData") # object: empSnpsByChr <-
# chromosomes <- paste0("chr", 1:19)
# snpDensityData <- list()
# for(c in chromosomes) {
#     snpDensityData[[c]] <- snp.density(chr=c, empRows=TRUE, empPositions = empSnpsByChr[[c]])
# }
# names(snpDensityData)[1:19] <- paste0("chr", 1:19)

### there's no need to run this function unless we have new genotype data.
### until then, you can just load this R object, snpDensityData, for plotting.
# save(snpDensityData, file="/group/palmer-lab/AIL/qtlmapping/snpDensity.allData.RData")