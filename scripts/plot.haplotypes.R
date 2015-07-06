# testing out DOQTL genome plot code.

# one way to visualize missing genotypes is to fill in genotyping gaps with NAs
# and map a different state.col to NA
# all you would need here is a pair of SNPs flanking the gaps between typed SNPs
# to delineate a region of missing/untyped variation.

##### LOAD DATA FROM PALMER PC ------------------------------------------------

setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing")
# chr19 haplotypes expressed as Lg or Sm ancestry (LL=2, LS=1, SS=0)
load("./dataFiles/chr19.hap.RData")
# genotyped relatives
load("./pedigree/genotypedRelatives.RData")
# genotyped sample ids
ids <- read.table("./genotyped.samples.txt", as.is=T)$V1
# autosome lengths in Mb
chrlen <- read.table("./dataFiles/chrLengths_mm10.txt",as.is=T)$V2 / 1e6

# snp positions
positions           <- read.table("./dataFiles/chr19.txt.legend.txt", header=F)[1:2]
names(positions)    <- c("snp", "pos")
typed.snps          <- read.table("./dataFiles/chr19.filtered.dosage", header=F, as.is=T)[1]
names(typed.snps)   <- c("snp")
positions           <- positions[positions$snp %in% typed.snps$snp,]
positions           <- positions$pos/1e6
rm(typed.snps)

# # distance between snps - NOT USED... FOR NOW.
# dist <- c()
# for (i in seq_along(positions[1:7882])){
#     dist[[i]] <- positions[i+1] - positions[i]
# }
# # currently no SNPs with > 1.85 Mb distance between them on chr19. This may
# # change once we filter out poorly imputed SNPs,excess missingness, etc. it will
# # also differ for chrs 10, 6, and 14, which contain regions where no SNPs were
# # typed by GBS.
# gaps<- which(dist > 5)
#
#
#
# ##### PLOT DATA FROM CHR 19 ---------------------------------------------------
#
# # identifiers for each haplotype state
# colors <- cbind(state = c("LL", "LS", "SS"),
#                 id = c("LG", "LGSM", "SM"),
#                 color = c("dodgerblue1", "goldenrod1", "firebrick1"))
#
# colors <- cbind(state = c("LL", "LS", "SS"),
#                 id = c("LG", "LGSM", "SM"),
#                 color = c("dodgerblue1", "goldenrod1", "firebrick1"))
#
# # match state colors to individual haplotype in chr19.hap
# state.cols <- colors[match(chr19.hap[[100]], colors[,1]),3]
#
# # set plot parameters
# par(plt = c(0.08, 0.95, 0.05, 0.88))
# chr.skeletons(chr=1:19, chrlen=chrlen[1:19])
#
# breaks <- c(1, which(state.cols[1:(length(state.cols)-1)] != state.cols[2:length(state.cols)]) + 1, length(state.cols))
#
# st.cols <- state.cols[breaks]
#
# # coordinates
# c <- 19
# xl <- c - 0.25
# xr <- c + 0.25
# y <- positions
#
# # draw the first haplotype
# rect(xl, y[breaks[1]], xr, y[breaks[2]],
#      col = st.cols[1], density=NA, border = st.cols[1])
#
# # draw subsequent haplotypes
# for (i in 2:length(st.cols)) {
#     rect(xl, y[breaks[i]], xr, y[breaks[i+1]],
#          col = st.cols[i], density=NA, border = st.cols[i])
# }
#
# legend(x = 17, y = 190, legend = colors[,2], col = colors[,3], pch = 15,
#        bg = "white")


############## TEST ------------------------------------------------------------

plot.LgSm.haplos <- function ( chr = 19, mice = ,
                               twoColors = c("dodgerblue1", "firebrick1"),
                               region = NULL,
                          #genoSamples = "./genotyped.samples.txt",
                          hapStates = "./dataFiles/chr19.hap.RData",
                          chrLengths = "./dataFiles/chrLengths_mm10.txt",
                          genoFile = "./dataFiles/chr19.filtered.dosage",
                          legend = "./dataFiles/chr19.txt.legend.txt") {

    # load haplotype states
    load(hapStates)
    # get chromosome lengths
    chrLens <- read.table(chrLens,as.is=T)$V2 / 1e6
    chrlength <- chrLens[chr]

    # get Mb positions of genotyped SNPs in the haplo data set
    positions <- read.table(legend, header=F)[1:2]
    typed.snps <- read.table(genoFile, header=F, as.is=T)[1]
    positions <- positions[positions$V1 %in% typed.snp$V1,]
    positions <- positions$pos/1e6
    rm(typed.snps)

    # to get 2 chromosomes for each mouse, i have to copy the vector of haps.
    # in each case, all 2 and 0 states (LG and SM, respectively) stay the same.
    # in one copy, all 1 states change to 2 and in the other copy, all 1 states
    # change to 0 to represent heterozygous genotypes on the 2 chromosomes.
    chr19.hap <- lapply(chr19.hap, as.data.frame)
    chr19.tmp <- chr19.hap
    #hap.tmp <- lapply(paste0("chr", chr, ".hap"), as.data.frame)
    hets <- c()
    for (i in seq_along(chr19.hap)) {
        chr19.hap[[i]][2] <- chr19.hap[[i]][1]
        names(chr19.hap[[i]]) <- c("one", "two")
        hets[[i]] <- which(chr19.hap[[i]][1] == 1)
        chr19.hap[[i]][[1]] <- replace(chr19.hap[[i]][[1]], hets[[i]],
                                     values=rep(2, times=length(hets[[i]])))
        chr19.hap[[i]][[2]] <- replace(chr19.hap[[i]][[1]], hets[[i]],
                                     values=rep(0, times=length(hets[[i]])))
        }

    # assign colors to LG(2) and SM(0)
    colors <- cbind.data.frame(state = c(2,0),
                    id = c("LG/J", "SM/J"),
                    color = twoColors)
    colors$color <- as.character(colors$color)

#     # get data for sib pairs
#     sib.haps <- c()
#     for (i in seq_along(sibPairs)) {
#         # prepare sib data
#         sib.haps[[i]] <- chr19.hap[names(chr19.hap) %in% sibPairs[[i]] ]
#     }
#
#     for (i in seq_along(sib.haps)){
#         for (mouse in seq_along(sib.haps[[i]])) {
#             sib.haps[[i]][[mouse]]$lcol <- colors[match(sib.haps[[i]][[mouse]]$one,
#                                                         colors[,1]),3]
#             sib.haps[[i]][[mouse]]$rcol <- colors[match(sib.haps[[i]][[mouse]]$two,
#                                                         colors[,1]),3]
#         } # for mouse
#     } # for i
#
#
#     # get breaks for left and right chromosomes
#     lbreaks <- vector("list", length=length(sib.haps))
#     rbreaks <- vector("list", length=length(sib.haps))
#     for (el in seq_along(lbreaks)){
#         lbreaks[[el]] <- list()
#         rbreaks[[el]] <- list()
#     }
#
#     for (i in seq_along(lbreaks)){
#         for (mouse in seq_along(sib.haps[[i]])) {
#             lbreaks[[i]][[mouse]] <- c(1, which(sib.haps[[i]][[mouse]]$lcol[1:length(sib.haps[[i]][[mouse]]$lcol)-1] != sib.haps[[i]][[mouse]]$lcol[2:length(sib.haps[[i]][[mouse]]$lcol)]) + 1, length(sib.haps[[i]][[mouse]]$lcol))
#             rbreaks[[i]][[mouse]] <- c(1, which(sib.haps[[i]][[mouse]]$rcol[1:length(sib.haps[[i]][[mouse]]$rcol)-1] != sib.haps[[i]][[mouse]]$rcol[2:length(sib.haps[[i]][[mouse]]$rcol)]) + 1, length(sib.haps[[i]][[mouse]]$rcol))
#         }
#     }
#
#     # map break points to color
#     sib.tmp <- vector("list", length=length(sib.haps))
#
#     for (i in seq_along(sib.haps)) {
#         sib.tmp[[i]] <- vector("list", length=length(sib.haps[[i]]))
#         for (j in seq_along(sib.haps[[i]])) {
#             sib.tmp[[i]][[j]] <- vector("list")
#             sib.tmp[[i]][[j]][["lcol"]] <- sib.haps[[i]][[j]]$lcol[lbreaks[[i]][[j]] ]
#             sib.tmp[[i]][[j]][["rcol"]] <- sib.haps[[i]][[j]]$rcol[rbreaks[[i]][[j]] ]
#         }
#         names(sib.tmp[[i]]) <- names(sib.haps[[i]])
#     }



        # set up plot area
        par(font = 2, font.axis = 2, font.lab = 2, las = 1,
            plt = c(0.2, 0.9, 0.12, 0.95))
        plot(1, 1, col = 0, xlim = c(1, (2 * length(sib.haps[[1]]) + 1)),
             ylim = c(0, chrlength+5), xaxt = "n", yaxs = "i", xlab = "",
             ylab = paste0("Chr ", chr, " position (Mb)"))

        for (i in 1:length(sib.haps[[1]])) {
            lines(2 * c(i,i), c(0, chrlength))
        }

        abline(h = 0:20 * 5, col = "grey80")
        #mtext(side = 1, at = 1.5, text = paste(names(sib.haps[[1]])))
        offset = 0.6
        for (c in 1:length(sib.haps[[1]])) {
            xm = 2 * c
            xl = xm - offset
            xr = xm + offset
            y = positions


#             # Plot sib2 (right) chromosome.
#             rect(1, y[sib.haps, xr, y[rbreaks[2]],
#                  col = rcol[1], density = NA, border = rcol[1])
        }


#         # coordinates
#         c <- 1
#         xl <- c - 0.11
#         xr <- c + 0.11
#         y <- positions




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
# matching to position


###### VERSION 2: HAP TO COLOR ------------------------------- 7/6/15 nmg
hap.to.color <- function (chr) {
    # read in hap file w/ user-defined path
    # hapFile <- paste0("chr", chr, ".hap")

    chr19.hap <- lapply(chr19.hap, as.data.frame)

    hets <- c()
    for (i in seq_along(chr19.hap)) {
        chr19.hap[[i]][2] <- chr19.hap[[i]][1]
        names(chr19.hap[[i]]) <- c("one", "two")
        hets[[i]] <- which(chr19.hap[[i]][1] == 1)
        chr19.hap[[i]][[1]] <- replace(chr19.hap[[i]][[1]], hets[[i]],
                                       values=rep(2, times=length(hets[[i]])))
        chr19.hap[[i]][[2]] <- replace(chr19.hap[[i]][[1]], hets[[i]],
                                       values=rep(0, times=length(hets[[i]])))
    }
    # assign colors to LG(2) and SM(0)
    #         colors <- cbind.data.frame(state = c(2,0),
    #                                    id = c("LG/J", "SM/J"),
    #                                    color = twoColors)
    #         colors$color <- as.character(colors$color)

    for (i in seq_along(chr19.hap)){
        for (chrom in seq_along(chr19.hap[[i]])) {
            lg <- which(chr19.hap[[i]][[chrom]] == 2)
            sm <- which(chr19.hap[[i]][[chrom]] == 0)
            chr19.hap[[i]][[chrom]] <- replace(chr19.hap[[i]][[chrom]], lg,
                                               rep("dodgerblue1", times=length(lg)))
            chr19.hap[[i]][[chrom]] <- replace(chr19.hap[[i]][[chrom]], sm,
                                               rep("firebrick1", times=length(sm)))
        } # for chrom
    } # for i

    # GET BREAK POINTS FOR LEFT AND RIGHT CHROMOSOMES
    breaks <- vector("list", length=2)
    names(breaks) <- c("left", "right")
    for (mouse in seq_along(breaks)) {
        breaks[[mouse]] <- vector("list", length=length(chr19.hap))
    }

    for (i in seq_along(breaks[["left"]])){
        breaks[["left"]][[i]] <- c(1, which(chr19.hap[[i]]$one[1:length(chr19.hap[[i]]$one)-1] != chr19.hap[[i]]$one[2:length(chr19.hap[[i]]$one)])+1, length(chr19.hap[[i]]$one))
    }
    for (i in seq_along(breaks[["right"]])){
        breaks[["right"]][[i]] <- c(1, which(chr19.hap[[i]]$two[1:length(chr19.hap[[i]]$two)-1] != chr19.hap[[i]]$two[2:length(chr19.hap[[i]]$two)])+1, length(chr19.hap[[i]]$two))
    }
    names(breaks[["left"]]) <- names(chr19.hap)
    names(breaks[["right"]]) <- names(chr19.hap)


    ### MAP BREAK POINTS TO COLORS
    col.breaks <- vector("list", length=2)
    names(col.breaks) <- c("left", "right")
    for (mouse in seq_along(col.breaks)) {
        col.breaks[[mouse]] <- vector("list", length=length(chr19.hap))
    }
    for (i in seq_along(chr19.hap)) {
        col.breaks[["left"]][[i]] <- chr19.hap[[i]]$one[breaks$left[[i]] ]
        col.breaks[["right"]][[i]] <- chr19.hap[[i]]$two[breaks$right[[i]] ]
    }
    names(col.breaks[["left"]]) <- names(chr19.hap)
    names(col.breaks[["right"]]) <- names(chr19.hap)


#------------------ POTENTIALLY THE END OF THIS FUNCTION.  RETURN COL.BREAKS


    # get data for sib pairs: test data set (family 1 in sib pairs)
    sibtest <- c(col.breaks$left[names(col.breaks$left) %in% sibPairs[[1]]],
                 col.breaks$right[names(col.breaks$right) %in% sibPairs[[1]]],
                 breaks$left[names(breaks$left) %in% sibPairs[[1]]],
                 breaks$right[names(breaks$right) %in% sibPairs[[1]]] )


    par(font = 2, font.axis = 2, font.lab = 2, las = 1,
        plt = c(0.2, 0.9, 0.12, 0.95))
    plot(1, 1, col = 0, xlim = c(1, (2 * length(sibtest) + 1)),
         ylim = c(0, chrlength+2), xaxt = "n", yaxs = "i", xlab = "",
         ylab = paste0("Chr ", chr, " position (Mb)"))

    for (a in 1:length(sibtest)) {
        lines(2 * c(a, a), c(0, chrlength))
    }

    abline(h = 0:20 * 5, col = "grey80")
    offset = 0.6

    for (c in 1:length(sibtest)) {
        xm = 2 * 1
        xl = xm - offset
        xr = xm + offset
        y = positions


        # need to structure the data to run a nested for loop
        # Plot sib2 (right) chromosome.
        rect(xl, y[breaks$left[[1]][1] ], xm, y[breaks$left[[1]][2] ],
             col = sibtest[[1]]$one[1], density = NA, border = sibtest[[1]]$one[1])

        rect(xm, y[breaks$right[[1]][1] ], xr, y[breaks$right[[1]][2] ],
             col = sibtest[[1]]$two[1], density = NA, border = sibtest[[1]]$two[1])

        rect(xl, y[breaks$left[[1]][i] ], xm, y[breaks$left[[1]][i+1] ],
             col = sibtest[[1]]$one[i], density = NA, border = sibtest[[1]]$one[i])
        rect(xm, y[breaks$right[[1]][i] ], xr, y[breaks$right[[1]][i+1] ],
             col = sibtest[[1]]$two[i], density = NA, border = sibtest[[1]]$two[i])
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



