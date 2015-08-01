# ADAPTED FROM DOQTL/EXTRACT.RAW.DATA - AUTHOR: DAN GATTI
################################################################################
# Extract the intensity or genotype data from the raw GeneSeek data files.
# Arguments: in.path: character vector of full paths to the GeneSeek data
#                     directories.
#            prefix: character vector of same length as in.path containing
#                    a prefix to add to each sample ID in data sets being
#                    processed.
#            out.path: character, full path to output directory.
#            simple.geno: if false, get genotypes in AT, AA, TG... format. if
#                    true, use the original DOQTL format (e.g. GG for G and H for
#                    any heterozygote class)
# Output: Writes out call.rate.batch.txt and genotypes.txt for all options.
#         Writes out x.txt and y.txt when type = "intensity".
################################################################################
extract.raw.data = function(in.path = ".", prefix, out.path = ".",
                            simple.geno = FALSE) {

    if(!missing(prefix)) {
        if(length(in.path) != length(prefix)) {
            stop(paste("extract.raw.data: The 'in.path' and 'prefix' matrices must",
                       "be the same length. Each prefix will be added to the matching",
                       "data set in in.path."))
        } # if(length(in.path) != length(prefix))
    } # if(!missing(prefix))
    tmp = sub("/$", "", in.path)
    if(!all(file.exists(tmp))) {
        stop(paste("extract.raw.data: Some of the in.path directories do not exist.",
                   paste(in.path[!file.exists(tmp)], collapse = ",")))
    } # if(!all(file.exists(in.path)))

    # Load gigasnp name/info file from JAX and remove the first 187 SNPs. These
    # are control sequences and weren't returned in our data. This should leave
    # you with a list of X snps, where X is the value in the 'Num Snps' field in
    # the header of the manifest (143,259 as opposed to Total Snps = 143,446).
   # load("gigasnps.RData")

    load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
    snps <- GM_snps[188:length(gigasnps$Name),]
    names(snps)[1] <- "SNP_ID" # so as not to be confused with Name in sampfile

    # Write out headers for the files.  This will overwrite existing files.
    x_file = file(paste(out.path, "x.txt", sep = "/"), open = "w")
    y_file = file(paste(out.path, "y.txt", sep = "/"), open = "w")
    g_file = file(paste(out.path, "geno.txt", sep = "/"), open = "w")
    writeLines(text = snps$SNP_ID[-length(snps$SNP_ID)], con = x_file, sep = "\t")
    writeLines(text = snps$SNP_ID[length(snps$SNP_ID)],  con = x_file, sep = "\n")
    writeLines(text = snps$SNP_ID[-length(snps$SNP_ID)], con = y_file, sep = "\t")
    writeLines(text = snps$SNP_ID[length(snps$SNP_ID)],  con = y_file, sep = "\n")
    writeLines(text = snps$SNP_ID[-length(snps$SNP_ID)], con = g_file, sep = "\t")
    writeLines(text = snps$SNP_ID[length(snps$SNP_ID)],  con = g_file, sep = "\n")
    call.rate.batch = NULL

    for(i in 1:length(in.path)) {
        # Get the sample IDs from the Sample_Map.txt file.
        samplefile = dir(path = in.path[i], pattern = "Sample_Map.txt", full.names = TRUE)
        # If not found, then quit.
        if(length(samplefile) == 0) {
            stop(paste("No file called 'Sample_Map.txt' was found in directory",
                       in.path[i], ".  Please make sure that the Sample_Map file is unzipped and",
                       "in the specified directory."))
        } # if(length(samplefile) == 0)
        samples = read.delim(samplefile, stringsAsFactors = FALSE)$Name

        # Find a file with "FinalReport" in the filename.
        rawfile = dir(path = in.path[i], pattern = "FinalReport", full.names = TRUE)
        # rawfile = dir(path = in.path, pattern = "FinalReport", full.names = TRUE)
        rawfile = rawfile[grep("txt", rawfile)]
        # If not found, then quit.
        if(length(rawfile) == 0) {
            stop(paste("No file with 'FinalReport' in the filename was found in directory",
                       in.path[i], ".  Please make sure that the FinalReport file is unzipped and",
                       "in the specified directory."))
        } # if(length(rawfile) == 0)
        # If there is more than one FinalReport file, then quit.
        if(length(rawfile) > 1) {
            stop(paste("There is more than one file with FinalReport in the filename.",
                       "Please place only one data set in each directory."))
        } # if(length(rawfile) > 1)

        # Read in the first sample.  The current format requires us to skip 9 lines
        # because there is no comment delimiter at the top of the file.
        print(paste("Reading", rawfile, "..."))
        rawfile = file(rawfile, open = "r")
        data = readLines(con = rawfile, n = 10)
        hdr = strsplit(data, split = "\t")
        cn = hdr[[10]]
        # Verify that we have all of the column names that we expect.
        column.names = c("SNP Name", "Sample ID", "Allele1 - Forward",
                         "Allele2 - Forward", "Allele1 - Top", "Allele2 - Top",
                         "Allele1 - AB", "Allele2 - AB", "GC Score", "X", "Y")
        columns = match(column.names, cn)
        if(any(is.na(columns))) {
            stop(paste("All of the expected column names were not found in the",
                       "FinalReport file. The missing column(s) are:", paste(
                           colnames[is.na(columns)], collapse = ",")))
        } # if(any(is.na(columns)))

        cr = rep(0, length(samples))
        samples.in.data = rep("", length(samples))

        # We read the files in and write them out to conserve memory. As computers
        # get larger, we may be able to keep everything in memory.
        for(j in 1:length(samples)) {
            # Read in the data for one sample.
            data = readLines(con = rawfile, n = nrow(snps))
            data = strsplit(data, split = "\t")
            data = matrix(unlist(data), nrow = length(data[[1]]), ncol = nrow(snps))
            dimnames(data) = list(cn, data[1,])
            samples.in.data[j] = data[rownames(data) == "Sample ID",1]
            print(paste("Sample", j, "of", length(samples), ":", samples.in.data[j]))
            if(!missing(prefix)) {
                samples.in.data[j] = paste(prefix[i], samples.in.data[j], sep = "")
            } # if(!missing(prefix))

            # Sort the data to match the SNP order.
            data = data[,match(snps$SNP_ID, colnames(data))]
            # X
            writeLines(samples.in.data[j], con = x_file, sep = "\t")
            xint = data[rownames(data) == "X",]
            writeLines(xint[-length(xint)], con = x_file, sep = "\t")
            writeLines(xint[length(xint)], con = x_file)
            # Y
            writeLines(samples.in.data[j], con = y_file, sep = "\t")
            yint = data[rownames(data) == "Y",]
            writeLines(yint[-length(yint)], con = y_file, sep = "\t")
            writeLines(yint[length(yint)], con = y_file)

            # Genotype: if simple.geno = FALSE, use two-letter genotypes
            geno = paste(data[rownames(data) == "Allele1 - Forward",],
                         data[rownames(data) == "Allele2 - Forward",], sep = "")
            cr[j] = mean(geno != "--")

            # Otherwise get one-letter genotypes (original DOQTL version)
            if (is.true(simple.geno)) {
                alleles = unique(geno)
                hets = alleles[!alleles %in% c("--", "AA", "CC", "GG", "TT")]
                geno[geno == "AA"] = "A"
                geno[geno == "CC"] = "C"
                geno[geno == "GG"] = "G"
                geno[geno == "TT"] = "T"
                geno[geno == "--"] = "-"
                geno[geno %in% hets] = "H"
            } # if (is.true(simple.geno))
            writeLines(samples.in.data[j], con = g_file, sep = "\t")
            writeLines(geno[-length(geno)], con = g_file, sep = "\t")
            writeLines(geno[length(geno)], con = g_file)
        } # for(j)
        close(rawfile)
        names(cr) = samples.in.data
        call.rate.batch = rbind(call.rate.batch, cbind(sample = names(cr),
                                                       call.rate = cr, batch = in.path[i]))
    } # for(i)
    close(x_file)
    close(y_file)
    close(g_file)
    colnames(call.rate.batch) = c("sample", "call.rate", "batch")
    write.table(call.rate.batch, paste(out.path, "call.rate.batch.txt", sep = "/"),
                sep = "\t", row.names = FALSE)
} # extract.raw.data()

extract.raw.data(in.path=".", out.path=".")

###############

# get intersection of gbs and giga SNP sets
# here i am using only empirically called GBS SNPs (within-sample imputation only)
# which were called by Shyam using UnifiedGenotyper
gbs.snps <- read.table("/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chrAll.filtered.snpinfo", sep="\t", header=F)[2:3]

# just get chr, pos, and name from the gigaMuga SNP map. I'm not sure if they're
# using 1-base or 0-base positions. My guess is one-base but I'm adding a 0-base
# column just in case. I'm reordering columns to match BED file format so I can
# feed the data to LiftOver (UCSC Genome Browser map version converter) in Galaxy.
# gigaMUGA SNPs are in mm9 format and ours are mm10.
giga.snps <- read.table("./SNP_Map.txt", sep="\t", header=T)[c(1,3,4)]
giga.snps$start <- giga.snps$Position -1
giga.snps$Chromosome <- paste0("chr", giga.snps$Chromosome)
giga.snps <- giga.snps[c(2,4,3,1)]
names(giga.snps)[3] <- "end"
write.table(giga.snps, file="../../LgSm-DataProcessing/dataFiles/giga.mm9.bed", sep="\t", row.names=F, col.names=T, quote=F)

# load the list of giga SNPs lifted over from mm9 to mm10
gmm10 <- read.table(file="./giga.mm9to10.mapped.bed", sep="\t", as.is=T)
names(gmm10) <- c("chr", "start", "end", "index")
gmm10$bedName <- paste0(gmm10$chr, ":", gmm10$end)

# get subset of giga snps that are represented in gbs data
validationSnps <- read.table("./SNP_Map.txt", sep="\t", header=T)
names(validationSnps)[1] <- "index"
val.tmp <- merge(validationSnps, gmm10, by="index", all.y=TRUE) # 142969 x 13
val.tmp$start <- NULL
val.tmp <- val.tmp[!duplicated(val.tmp$bedName),] #140268

# gbs snps that are found in giga snp list
gbs.in.giga <- which(gbs$end %in% val.tmp$end)
gbs <- gbs[gbs.in.giga,]
names(gbs)[4] <- "bedName"
gbs$start <- NULL

# SNP_map for gigaSNPs found in gbs with both mm9 and mm10 positions listed
testVal <- val.tmp[which(val.tmp$end %in% gbs$end),]
names(testVal)[4] <- "mm9.pos"
names(testVal)[11] <- "mm10.pos"
testVal$chr <- NULL
validationSet <- testVal

# match giga samples to row indexes in genotyped.samples.txt so I can pull out
# their gbs genotypes later. they'll still need to be matched to the actual
# bases and Lg/Sm
gigasamp <- read.table("Sample_Map.txt", sep="\t", header=T)$Name
allsamp<- read.table("../../LgSm-DataProcessing/genotyped.samples.txt", sep="\t", header=F)$V1
sampidx <- which(allsamp %in% gigasamp)

# snp list to search filtered.dosage files
snplist <- paste0("ail.chr", validationSet$Chromosome,".", validationSet$mm10.pos)

save.image()
### next use haplo scripts to read in genotype info for gigasamples at validationSet snps


##################### PART 2 ################################
# get nucleotide and Lg/Sm ancestry information on giga snps
is.lgsm.giga <- function(chromosome,
                         legendPathPrefix = "/group/palmer-lab/AIL/knownSNPs/imputeHaplotypes/",
                         haploPathPrefix = "/group/palmer-lab/AIL/knownSNPs/imputeHaplotypes/",
                         filtDosPathPrefix = "/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/")  {


    legendFile      <- paste0(legendPathPrefix, chromosome, ".txt")
    hapFile         <- paste0(haploPathPrefix, chromosome, ".hap")

    # FIND WHICH STRAIN THE ALTERNATIVE ALLELE IS FROM
    legend          <- read.table(legendFile, header=F, as.is=T)[c(1,3,4)] # anme, ref, alt
    haps            <- read.table(hapFile, header=F, as.is=T)
    legend          <- cbind(legend, haps)
    names(legend)   <-  c("name","ref","alt", "LG", "SM")
    rm(legendFile, hapFile, haps)

    load("/group/palmer-lab/AIL/knownSNPs/gigaMUGA/snplist.gigamice.RData") # obj = snplist

    legend <- legend[legend$name %in% snplist,]

    return(legend)

}
# run is.lgsm.giga
chr = paste0("chr", 1:19)
gbsingiga <- lapply(chr, is.lgsm.giga)
save(gbsingiga, "gbsingiga.RData")

# how many snps were giga genotyped per chromosome?
gigasnps.perChr <- c()
for (i in seq_along(gbsingiga)) {
    gigasnps.perChr <- c(gigasnps.perChr, length(gbsingiga[[i]]$name))
}
names(gigasnps.perChr) <- chr

# chr19 result
$ :'data.frame':       52 obs. of  5 variables:
    ..$ name: chr [1:52] "ail.chr19.5650883" "ail.chr19.9876922"...
..$ ref : chr [1:52] "T" "G" "A" "T" ...
..$ alt : chr [1:52] "C" "T" "G" "C" ...
..$ LG  : int [1:52] 0 1 1 0 1 1 1 0 0 0 ...
..$ SM  : int [1:52] 1 0 0 1 0 0 0 1 1 1 ...

# -----------------------------------------------------------------------------
# write giga snp names out to text files
chr <- paste0("chr", 1:19)
cmds <- c()
for (chr in seq_along(gbsingiga)) {
    write.table(as.data.frame(gbsingiga[[chr]]$name), file=paste0("chr", chr, ".giganames"),
                col.names=F, row.names=F, quote=F, sep="\t")
    cmds <- c(cmds, paste0("find /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr", chr,
                           ".filtered.dosage -type f -exec grep -iFf chr", chr,
                           ".giganames {} + > chr", chr, ".snps"))
}
write.table(cmds, file="find.gigasnps.cmds", sep="\t", col.names=F, row.names=F, quote=F)

# in unix shell, find only the desired SNPs and write them to new files
# find /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr19.filtered.dosage -type f -exec grep -iFf chr19.giganames {} + > chr19.snps

# now match the columns to individuals. there are 3 columns in each filtered.dosage
# file (snp name, ref allele, alt allele) so we must add 3 to the sample indexes in
# order to extract the right columns.
sampidx$idxPlus3 <- sampidx$idx + 3

gbsGenos <- c()
for (i in 1:19){
    filename <- paste0("./gigaGenos/chr", i, ".snps")
    gbsGenos[[i]] <- read.table(filename, sep=" ", header=F, as.is=T)[sampidx$idxPlus3]
}
names(gbsGenos) <- paste0("chr", 1:19)
for (i in seq_along(gbsGenos)) {
    names(gbsGenos[[i]]) <- sampidx$id
}

# ------------------------------------------------------------------------------

#  genos contains raw dosage genotypes for the relevant samples/snps; gbsGenos
# legend is SNP info table, e.g. gbsingiga
load("/group/palmer-lab/AIL/knownSNPs/gigaMUGA/gbsingiga.RData")
load("/group/palmer-lab/AIL/knownSNPs/gigaMUGA/snplist.gigamice.RData")

dosGiga <- lapply(chr=chr, convert.hap.to.allele)

convert.hap.to.allele <- function (chr, legend = gbsingiga,
                                   genos = gbsGenos,
                                   samples=sampidx$id
) {
    # load haplotype data and pull out gigaMice
    haploFile = paste0("/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers/",
                       chr, ".hap.RData")
    load(haploFile)

    genos <- get(paste0(chr, ".hap"))
    # rm(paste0(chr, ".hap"))

    #### IMPORTANT - CHANGE gbsGenos to GENOS
    # gbs dosage for a het is 0.5-1.5, the values Shyam used for CFWs
    genoClass <- c()
    for (mouse in 1:24){
        genoClass[[chr]][[mouse]] <- cut(genos[[chr]][[mouse]],
                                         breaks=c(0, 0.5, 1.5, 2),
                                         labels=c("R", "H", "A"), dig.lab=4,
                                         right=TRUE, include.lowest=TRUE)
        # snps in rows; ids in columns
    }
    genoClass[[chr]] <- lapply(genoClass[[chr]], as.character)
    genoClass[[chr]] <- do.call(cbind, genoClass[[chr]])
    colnames(genoClass[[chr]]) <- samples
    row.names(genoClass[[chr]]) <- legend[[chr]]$name

    # get nucleotides
    for (mouse in seq_along(genoClass[[chr]])) {
        for (i in seq_along(genoClass[[chr]][[mouse]])){
            if(genoClass[[chr]][[mouse]][i] == "A"){
                genoClass[[chr]][[mouse]][i] = legend[[chr]][i,"alt"]
            }
        }
    }
    for (mouse in seq_along(genoClass[[chr]])) {
        for (i in seq_along(genoClass[[chr]][[mouse]])){
            if(genoClass[[chr]][[mouse]][i] == "R"){
                genoClass[[chr]][[mouse]][i] = legend[[chr]][i,"ref"]
            }
        }
    }
    return(genoClass)
}





