# testing out DOQTL genome plot code.

# one way to visualize missing genotypes is to fill in genotyping gaps with NAs
# and map a different state.col to NA
# all you would need here is a pair of SNPs flanking the gaps between typed SNPs
# to delineate a region of missing/untyped variation.

### 1. Load data from Palmer PC ------------------------------------------------
setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing")

# chr19 haplotypes expressed as Lg or Sm ancestry
load("./dataFiles/chr19.hap.RData")

# autosome lengths in Mb
chrlen <- read.table("./dataFiles/chrLengths_mm10.txt",as.is=T)$V2 / 1e6
names(chrlen) <- 1:19

# genotyped samples
ids <- read.table("./genotyped.samples.txt", as.is=T)$V1
names(chr19.hap) <- ids

# snp positions
positions <- read.table("./dataFiles/chr19.txt.legend.txt", header=F)[1:2]
names(positions) <- c("snp", "pos")
typed.snps <- read.table("./dataFiles/chr19.filtered.dosage", header=F, as.is=T)[1]
names(typed.snps) <- c("snp")
positions <- positions[positions$snp %in% typed.snps$snp,]
positions <- positions$pos/1e6
rm(typed.snps)

# distance between snps
dist <- c()
for (i in seq_along(positions[1:7882])){
    dist[[i]] <- positions[i+1] - positions[i]
}
# currently no SNPs with > 1.85 Mb distance between them on chr19. This may
# change once we filter out poorly imputed SNPs,excess missingness, etc. it will
# also differ for chrs 10, 6, and 14, which contain regions where no SNPs were
# typed by GBS.
gaps<- which(dist > 5)

### 2. Get IDs of interest from ped --------------------------------------------

load("./pedigree/info.phenoSamples.RData")
info <- info[order(info$gen, info$dam),]
# get info for genotyped samples only (ids)
info$id <- info$id - 0.1
info <- info[info$id %in% ids,]

# retain only mice with siblings
sibs.tmp <- c()
sibs <- list()
for (i in 1:(length(info$dam)-1)) {
    if (info$dam[i] == info$dam[i+1]) {
        sibs.tmp <- rbind(sibs.tmp, info[i,], info[i+1,])
    } # gets df containing sib pair info
    if (info$dam[i] == info$dam[i+1]) {
        sibs[[i]] <- c(info[i,1], info[i+1,1])
    } # makes list of sibs for each dam
}
sibs <- sibs[!sapply(sibs, is.null)]
names(sibs) <- unique(sibs.tmp$dam)


### 3. Plot data ---------------------------------------------------------------

# identifiers for each haplotype state
colors <- cbind(state = c("LL", "LS", "SS"),
                id = c("LG", "LGSM", "SM"),
                color = c("dodgerblue1", "goldenrod1", "firebrick1"))

# match state colors to individual haplotype in chr19.hap
state.cols <- colors[match(chr19.hap[[100]], colors[,1]),3]

# set plot parameters
par(plt = c(0.08, 0.95, 0.05, 0.88))
chr.skeletons(chr=1:19, chrlen=chrlen[1:19])

breaks <- c(1, which(state.cols[1:(length(state.cols)-1)] != state.cols[2:length(state.cols)]) + 1, length(state.cols))

st.cols <- state.cols[breaks]

# coordinates
c <- 19
xl <- c - 0.25
xr <- c + 0.25
y <- positions

# draw the first haplotype
rect(xl, y[breaks[1]], xr, y[breaks[2]],
     col = st.cols[1], density=NA, border = st.cols[1])

# draw subsequent haplotypes
for (i in 2:length(st.cols)) {
    rect(xl, y[breaks[i]], xr, y[breaks[i+1]],
         col = st.cols[i], density=NA, border = st.cols[i])
}

legend(x = 17, y = 190, legend = colors[,2], col = colors[,3], pch = 15,
       bg = "white")


############## TEST


plot.haplos <- function ( chr = 19, mice = ,
                          #genoSamples = "./genotyped.samples.txt",
                          hapStates = "./dataFiles/chr19.hap.RData",
                          chrLens = "./dataFiles/chrLengths_mm10.txt",
                          genoFile = "./dataFiles/chr19.filtered.dosage",
                          legend = "./dataFiles/chr19.txt.legend.txt") {

#   get list of genotyped samples
#     all.mice <- read.table(genoSamples, as.is=T)$V1
#     names(all.mice) <- ids

    # load haplotype states
    load(hapStates)

    # get chromosome lengths
    chrlen <- read.table(chrLens,as.is=T)$V2 / 1e6

    # get Mb positions of genotyped SNPs in the haplo data set
    positions <- read.table(legend, header=F)[1:2]
    typed.snps <- read.table(genoFile, header=F, as.is=T)[1]
    positions <- positions[positions$V1 %in% typed.snp$V1,]
    positions <- positions$pos/1e6
    rm(typed.snps)

    # identifiers for each haplotype state
    colors <- cbind(state = c("LL", "LS", "SS"),
                    id = c("LG", "LGSM", "SM"),
                    color = c("dodgerblue1", "goldenrod1", "firebrick1"))
    colors$color <- as.character(colors$color)

    # plot

    for (dam in seq_along(sibs)) {

        # prepare sib data
        sibs.col <- data.frame(chr19.hap[which(names(chr19.hap) == sibs[[dam]][1])],
                               chr19.hap[which(names(chr19.hap) == sibs[[dam]][2])])
        names(sibs.col) <- c(sibs[[dam]][1], sibs[[dam]][2])
        lcol <- colors[match(sibs.col[[1]], colors[,1]),3]
        rcol <- colors[match(sibs.col[[2]], colors[,1]),3]
        lbreaks = c(1, which(lcol[1:(length(lcol)-1)] != lcol[2:length(lcol)]) +
                        1, length(lcol))
        rbreaks = c(1, which(rcol[1:(length(rcol)-1)] != rcol[2:length(rcol)]) +
                        1, length(rcol))
        lcol = lcol[lbreaks]
        rcol = rcol[rbreaks]

        # set up plot area
        par(font = 2, font.axis = 2, font.lab = 2, las = 1,
            plt = c(0.2, 0.9, 0.12, 0.95))
        plot(1, 1, col = 0, #xlim = c(1, (2 * length(chrlen[19]) + 1)),
             ylim = c(0, 65), xaxt = "n", yaxs = "i", xlab = "",
             ylab = "Mb position on chr19")
        lines(c(1,1), c(0,chrlen[19]))
        abline(h = 0:20 * 5, col = "grey80")
        mtext(side = 1, at = 1:length(chrlen), text = paste(names(sibs.col)[1], "  ",
                                                            names(sibs.col)[2]))

        # coordinates
        c <- 1
        xl <- c - 0.11
        xr <- c + 0.11
        y <- positions

        # Plot the sib1 (left) chromosome
        rect(xl, y[lbreaks[1]], 1, y[lbreaks[2]],
             col = lcol[1], density = NA, border = lcol[1])
        # Plot sib2 (right) chromosome.
        rect(1, y[rbreaks[1]], xr, y[rbreaks[2]],
             col = rcol[1], density = NA, border = rcol[1])


        for(i in lbreaks[2:length(lbreaks)]) {
            rect(xl, y[lbreaks[i]], 1, y[lbreaks[i+1]],
                 col = lcol[i], density = NA, border = lcol[i])
        }
        for (i in rbreaks[2:length(rbreaks)]){
            # Plot the right side of the chromosome.
            rect(1, y[rbreaks[i]], xr, y[rbreaks[i+1]],
                 col = rcol[i], density = NA, border = rcol[i])
        }


        legend(x = 1.24, y = 60, legend = colors[,2], col = colors[,3], pch = 15,
               bg = "white")

    }
}







########### Chromosome skeletons -----------------------------------------------
# Draw the skeleton of the chromosomes.
# Arguments: chr: vector with chr names to plot.
chr.skeletons = function(chr, chrlen, ...) {

    # Plot a coordinate system to draw on.
    par(font = 2, font.axis = 2, font.lab = 2, las = 1, plt = c(0.14, 0.95,
                                                                0.1, 0.9))
    plot(1, 1, col = 0, xlim = c(1, (length(chrlen) + 1)),
         ylim = c(0, 200), xaxt = "n", yaxs = "i", xlab = "",
         ylab = "Mb")

    # Draw chr skeletons
    for(i in 1:length(chrlen)) {
        lines(c(i, i),
              c(0, chrlen[i]))
    } # for(i)

    # Draw light lines at 10 Mb and heavier lines at 50 Mb.
    abline(h = 0:20 * 10, col = "grey80")
    abline(h = 0:4  * 50, col = "grey50", lwd = 1.2)
    mtext(side = 1, at = 1:length(chrlen), text = names(chrlen))

}

# GENOMIC POINTS - POTENTIALLY UNNECESSARY
genomic.points = function(chr = NULL, loc = NULL) {
    if(!is.null(chr) & !is.null(loc)) {
        if(length(chr) != length(loc)) {
            stop(paste("chr and loc are not the same length. Please use chr and",
                       "loc vectors of the same length."))
        } # if(length(chr) != length(loc))
        points(chr, loc, pch = "-", cex = 2, lwd = 2, col = 2)
    } # if(!is.null(chr) & !is.null(loc))
}
