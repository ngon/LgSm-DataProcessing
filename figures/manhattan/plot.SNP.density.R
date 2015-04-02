# plot SNP density by chr in ggplot2
# data frame should have SNPs in rows, and a column for chromosome, position,
# and the factor you want to group it by, if any (e.g. empirical vs. imputed).


chrLens    <- read.table("chrLengths_mm10.txt")$V2
chrLens    <- c(0, cumsum(chrLens*1.0))

labPos     <- read.table("chrLengths_mm10.txt")$V2/2.0
labPos     <- chrLens[1:19]+labPos

chromosomes <- paste0("chr", 1:19)



#########

# need a df with chr, position, is.empirical, is.imputed

plot.SNP.density <- function(chr, outfile=NULL,title=NULL) {

    if (is.null(outfile)) {
        outfile <- paste0(chr, ".SNPdensity.png")
    }
    if (is.null(title)) {
        title <- paste0("SNP density on ", chr)
    }
    png(file=outfile, height=400, width=500, bg="transparent")
    #pvals <- c()
    positions <- c()
    #cols <- c()

    for (chrom in 1:19) {
        filename <- paste0("chr", chrom, ".filtered.dosage.txt")
        print(paste("Starting chromosome ", chrom))
        output <- read.table(filename, header=T, as.is=TRUE)
        #output$chr <- as.numeric(sub("chr", "", output$chr))
        #pvals <- c(pvals, -log10(output$p_lrt))
        positions <- c(positions, chrLens[output$chr]+output$ps)
#         if (chrom %% 2) {
#             cols <- c(cols, rep(oddcolor, nrow(output)))
#         } else {
#             cols <- c(cols, rep(evencolor, nrow(output)))
#         }
        print(paste("Done with chromosome", chrom))
    }
    ##### plotting code goes here
}
#########


snps<-read.table("snp132_ucsc_hg19.bed",sep="\t",header=F)
colnames(snps)<-c("chr","start","end","id","score","strand")

# put the chromosomes in order
goodChrOrder <- paste("chr",c(1:19,"X","Y"),sep="")
snps$chr <- factor(snps$chr,levels=goodChrOrder)

# Plot the densities of snps in the bed file for each chr seperately
library(ggplot2)
snpDensity<-ggplot(snps) +
    geom_histogram(aes(x=start),binwidth=1e6) + # pick a binwidth that is not too small
    facet_wrap(~ chr,ncol=2) + # seperate plots for each chr, x-scales can differ from chr to chr
    ggtitle("Density of SNP-132 across hg19") +
    xlab("Position in the genome") +
    ylab("SNP density") +
    theme_bw() # I prefer the black and white theme

##### TO GROUP PLOTS BY A FACTOR, USE FACET_GRID
# ggplot(yourData) +
#   geom_histogram(aes(x=position),binwidth=1e6) +
#    facet_grid(group ~ chr)

# save the plot to .png file
png("snp132_density.png",500,750)
print(snpDensity)
dev.off()