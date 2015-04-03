#### plot.SNP.density ######
#### a function for plotting SNP density for each chromosome

# Run this in /group/palmer-lab/AIL/GBS/genoSummaries
# to do: read in each chromosome at once, get the information you need and then
# discard from memory.
# all you need is ps col from each chrInfo file. immediately divide all the ps
# by 1 Mb (or whatever bin size you want) and store the results AS INTEGERS.
# then: chr1bins<- tapply(sort(unique(chrInfoTable)), chrInfoTable, count)
# this will give you the number of SNPs in each 1 Mb bin. save this into a list
# use hist(chr1bins, plot=FALSE)
# then match by row name to list of empirical snps

# get preimpute filenames except for chrX (idx 20), which is empty
preimpFiles <- list.files(path="../preimpute/", pattern="*.preimpute.geno")
preimpFiles <- preimpFiles[-20]
# extract row names from each file
get.emp.rownames <- function(file) {read.table(file, sep="\t", header=F)[3]}
empSnps <- lapply(file.path("../preimpute/", preimpFiles), get.emp.rownames)


plot.SNP.density <- function(infoFile, empFile) {

    print("Finding the necessary files...")
    chromosomes <- paste0("chr", 1:19)
    filenames <- list()
    for (i in chromosomes){
        filenames[i] <- paste0("../dosage/", i, ".filtered.snpinfo")
    }

    read.info <- function(filename) {read.table(filename, sep="\t", header=F,
                                     col.names=c("snp", "ps", "chr"))}
    info <- lapply(filenames, read.info)
    info <- do.call(what=rbind.data.frame, info)

    # make new column marking each SNP as empirical=TRUE or FALSE
    emp<- read.table("./empirical.snp.list.txt", sep="\t", header=T)[1]
    info$is.emp <- info$snp %in% emp$snp

    p<- ggplot(data=info, aes(x=ps/1e6)) +
        geom_histogram(aes(x=ps/1e6), binwidth=1, fill="goldenrod2") +
        geom_histogram(data=info[info$is.emp == T,],
                       binwidth=1, fill="steelblue3") +
        facet_grid(chr~.) +
        ggtitle("Genome-wide SNP density") +
        xlab("Chromosome position (Mb)") +
        ylab("SNP density") +
        theme_bw() +
        theme(axis.title.x=element_text(size=12),
              axis.title.y=element_text(size=12),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10))

    save(p, "./snpDensity.RData")
    #     png(file=paste0(outfile, ".png"), height = 11, width=8, units="in")
    #     print(p)
    #     dev.off()

    }





