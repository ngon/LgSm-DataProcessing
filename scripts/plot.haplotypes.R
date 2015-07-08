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

# assign colors to LG(2) and SM(0)
# colors <- cbind.data.frame(state = c(2,0),
#                            id = c("LG/J", "SM/J"),
#                            color = c(lg.col, sm.col))
# colors$color <- as.character(colors$color)



######  HAP TO COLOR ------------------------------- 7/6/15 nmg
# Haplotype list must be preloaded before running hap.to.color. I have stored
# these as RData files for each chromosome. Name format is "chrN.hap.RData" and
# the object within is always "chrN.hap". It returns two lists: breaks, with
# list positions where a haplotype switch occurs, and col.breaks, where breaks
# is matched to colors representing LG and SM.
load("./dataFiles/chr19.hap.RData")
hapdata <- hap.to.color(hap=chr19.hap)

hap.to.color <- function (hap, lg.col="dodgerblue1", sm.col="firebrick1") {

    hap <- lapply(hap, as.data.frame)

    # make a left and right copy of each chromosome. for all heterozygous
    # states (1), change the left chr to LG(2) and the right chr to SM(0).
    hets <- c()
    for (i in seq_along(hap)) {
        hap[[i]][2] <- hap[[i]][1]
        names(hap[[i]]) <- c("one", "two")
        hets[[i]] <- which(hap[[i]][1] == 1)
        hap[[i]][[1]] <- replace(hap[[i]][[1]], hets[[i]],
                                       values=rep(2, times=length(hets[[i]])))
        hap[[i]][[2]] <- replace(hap[[i]][[1]], hets[[i]],
                                       values=rep(0, times=length(hets[[i]])))
    }
    rm(hets) # cleanup

    # now replace LG(2) and SM(0) states with lg.col and sm.col (color names)
    for (i in seq_along(hap)){
        for (chrom in seq_along(hap[[i]])) {
            lg <- which(hap[[i]][[chrom]] == 2)
            sm <- which(hap[[i]][[chrom]] == 0)
            hap[[i]][[chrom]] <- replace(hap[[i]][[chrom]], lg,
                                               rep(lg.col, times=length(lg)))
            hap[[i]][[chrom]] <- replace(hap[[i]][[chrom]], sm,
                                               rep(sm.col, times=length(sm)))
        } # for chrom
    } # for i

    #  set up breaks list structure
    breaks <- vector("list", length=2)
    names(breaks) <- c("left", "right")
    for (mouse in seq_along(breaks)) {
        breaks[[mouse]] <- vector("list", length=length(hap))
    }
    # find breakpoints for left and right chromosomes
    for (i in seq_along(breaks[["left"]])){
        breaks[["left"]][[i]] <- c(1,
        which(hap[[i]]$one[1:length(hap[[i]]$one)-1] != hap[[i]]$one[2:length(hap[[i]]$one)])+1,
        length(hap[[i]]$one))
    }
    for (i in seq_along(breaks[["right"]])){
        breaks[["right"]][[i]] <- c(1,
        which(hap[[i]]$two[1:length(hap[[i]]$two)-1] != hap[[i]]$two[2:length(hap[[i]]$two)])+1,
        length(hap[[i]]$two))
    }
    names(breaks[["left"]]) <- names(hap)
    names(breaks[["right"]]) <- names(hap)


    # set up col.breaks list structure
    col.breaks <- vector("list", length=2)
    names(col.breaks) <- c("left", "right")
    for (mouse in seq_along(col.breaks)) {
        col.breaks[[mouse]] <- vector("list", length=length(hap))
    }
    # map breakpoints to colors
    for (i in seq_along(chr19.hap)) {
        col.breaks[["left"]][[i]] <- hap[[i]]$one[breaks$left[[i]] ]
        col.breaks[["right"]][[i]] <- hap[[i]]$two[breaks$right[[i]] ]
    }
    names(col.breaks[["left"]]) <- names(hap)
    names(col.breaks[["right"]]) <- names(hap)

    list(breaks=breaks, col.breaks=col.breaks)
} # hap.to.color

# structure of the list is: listName [[breaks]] [["right"]] [[mouseID]]
# names(hapdata) <- c("breaks", "col.breaks")
# save memory by >  rm(chr19.hap)

### FUNCTION: plot.haplos ------------------------------------------------------
# You must have ids object loaded to run. This can be an element of sibPairs or
# fams, located in ./pedigree/genotypedRelatives.Rdata. Contents of each list
# are described in the README file in pedigree folder.
load("./pedigree/genotypedRelatives.RData")

plot.haplos <- function ( chr = 19, ids = sibPairs[[1]], breaks = hapdata[["breaks"]],
                          col.breaks = hapdata[["col.breaks"]], gen=NULL,
                               lg.col = "dodgerblue1", sm.col = "firebrick1",
                               region = NULL,
                               chrLengths = "./dataFiles/chrLengths_mm10.txt",
                               genoFile = "./dataFiles/chr19.filtered.dosage",
                               legend = "./dataFiles/chr19.txt.legend.txt",
                               outfilePrefix = "./figures/familyHaplos_",
                               pdfInches = 4) {

    # get chromosome length
    chrlength <- read.table(chrLengths,as.is=T)$V2[chr] / 1e6

    # get Mb positions of genotyped SNPs in the haplo data set
    # pos$V1 is SNP name to match to typed.snps$V1. pos$V2 = bp position.
    pos <- read.table(legend, header=F)[1:2]
    typed.snps <- read.table(genoFile, header=F, as.is=T)[1]
    pos <- pos[pos$V1 %in% typed.snps$V1,]
    pos <- pos$V2/1e6
    rm(typed.snps) # cleanup

    # load graphics device---------------------------
    pdf(file=paste0(outfilePrefix, gen, ".pdf"), width=pdfInches, height=pdfInches)

   # plot-------------------------------------------
    for (fam in seq_along(ids)) {
        # get data for sib pairs: test data set (family 1 in sib pairs)
        sibtest <- c(breaks$left[names(breaks$left) %in% ids[[fam]]],
                     breaks$right[names(breaks$right) %in% ids[[fam]]],
                     col.breaks$left[names(col.breaks$left) %in% ids[[fam]]],
                     col.breaks$right[names(col.breaks$right) %in% ids[[fam]]])

        # preserve mouse IDs
        attr(sibtest, "mouse_ids") <- unique(names(sibtest))
        # get half the length of sibset to divide the cols from breaks
        div = length(sibtest)/2

        # names for sibtest
        name.tmp <- paste0(names(sibtest), c(rep(".l", times=div/2), rep(".r", times=div/2),
                                             rep(".l", times=div/2), rep(".r", times=div/2)))
        name.bk <- paste0(name.tmp[1:div], rep(c("bk"), times=div))
        name.col <- paste0(name.tmp[div+1:div], rep(c("col"), times=div))
        names(sibtest) <- c(name.bk, name.col)
        rm(name.tmp, name.bk, name.col) # cleanup

        # set up plot area. use bold fonts and horizontal axis labels. these
        # parameters print well on a 4x4 inch graphics device.
        par(font = 2, font.axis = 2, font.lab = 2, las = 1,
            plt = c(0.2, 0.9, 0.12, 0.95), mar=c(2,4,1,2.7)+0.1)

        # draw chromosome skeletons and label the y-axis
        plot(1, 1, col = 0, xlim = c(1, (2 * (div/2) + 1)),
             ylim = c(0, chrlength+2), xaxt = "n", yaxs = "i", xlab = "", ### VAR chrlength
             ylab = paste0("Chr ", chr, " position (Mb)"), cex=0.9) ### VAR chr

        # draw chromosome skeletons
        for (a in 1:(div/2)) {
            lines(2 * c(a, a), c(0, chrlength))
        }
        # gridlines and a legend for clarity
        abline(h = 0:20 * 5, col = "grey80")

        # add a 'legend' (if you use the legend function in baseR, it takes up
        # a lot of space and the position may need to be adjusted for different
        # sizes of families. putting labels in the margin avoids this issue).
        mtext(side = 4, at = c((chrlength/2)+5, chrlength/2), adj = -0.2,
              text = c("LG/J", "SM/J"), col = c(lg.col, sm.col), cex = 0.8)

        # x-axis labels for each mouse in sibtest
        for (mouse in seq_along(attr(sibtest, "mouse_ids"))) {
            mtext(side = 1, at = mouse * 2, line = 0.2,
                  text = attr(sibtest, "mouse_ids")[mouse])
        }

        for (c in 1:(div/2)) { # for all elements containing LEFT breaks
            offset = 0.5
            xm = 2 * c
            xl = xm - offset
            xr = xm + offset

            for (idx in seq_along(sibtest[[c]])) {
                if (idx == 1){
                    rect(xl, pos[sibtest[[c]][1] ], xm, pos[sibtest[[c]][2] ],
                        col = sibtest[[c+div]][1], density = NA,
                         border = NA )
                    rect(xm, pos[sibtest[[(c+(div/2))]][1] ], xr, pos[sibtest[[c+(div/2)]][2] ],
                         col = sibtest[[(c+(div/2))+div]][1], density = NA,
                         border = NA )
                } else {
                    if (idx > 1) {
                        rect(xl, pos[sibtest[[c]][idx] ], xm, pos[sibtest[[c]][idx+1] ],
                             col = sibtest[[c+div]][idx], density = NA,
                             border = NA )
                        rect(xm, pos[sibtest[[(c+(div/2))]][idx] ], xr, pos[sibtest[[c+(div/2)]][idx+1] ],
                             col = sibtest[[(c+(div/2))+div]][idx], density = NA,
                             border = NA )
                    } # if idx > 1
                } # else
            } # for (idx in seq_along (sibtest[[c]]))
        } # for (c in 1:div/2)
    } # for (fam in seq_along(ids))
    dev.off()
} # plot.haplos

# i want to put the plots into several different PDF files to avoid making them
# too large to open. here i'm grouping them by generation. here i borrowed code
# from ./pedigree/genotyped.relatives.R: i load the pedigree and get info for
# mice that were genotyped and have siblings (e.g. mice in sibPairs). I make
# sure that the list of dams in info is identical to names(sibPairs). Now I can
# determine what range of sibPair elements comprises each generation of AILs
# and plot them in separate PDF files.
load("./pedigree/info.phenoSamples.RData")
info <- info[order(info$gen, info$dam),]
ids <- read.table("./genotyped.samples.txt", as.is=T)$V1
info$id <- info$id - 0.1
info <- info[info$id %in% ids,]

test <- info[info$dam %in% as.numeric(names(sibPairs)),c(1:3)]
test <- test[which(!duplicated(test$dam)),]
test <- test[order(test$dam),]
identical(test$dam, as.integer(names(sibPairs))) # yes

test$gen <- as.factor(test$gen)
row.names(test) <- 1:509
test <- split(test, test$gen, drop=T)
rm(info, ids, test) # cleanup

genrows <- lapply(test, row.names)
genrows <- lapply(genrows, as.integer)

chk <- sibPairs[seq_along(sibPairs) %in% genrows[["56"]]]

# first run: there are a few things i noticed.
# sometimes the chromosome arms appear to be different sizes. for example:
# plot #5, the right chromosome for 54380 is shorter.
# plot 10: R chr of 54867 is shorter
# plot 12: R chr of 54383 is shorter
# plot 26: R chr of 54830 is shorter

# missing ID for first mouse in plot 8 (sib of 54816)
# this is a quad - there are 4 sibs. the plot did not expand to include all
# of them. the label 54816 is the first mouse in the list.

plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["56"]]],
                      breaks = breaks, col.breaks = col.breaks, gen="56",
                          lg.col = "dodgerblue1", sm.col = "firebrick1",
                          region = NULL,
                          chrLengths = "./dataFiles/chrLengths_mm10.txt",
                          genoFile = "./dataFiles/chr19.filtered.dosage",
                          legend = "./dataFiles/chr19.txt.legend.txt",
                          outfilePrefix = "./figures/familyHaplos_",
                          pdfInches = 4)


chL <- breaks[["left"]][["54830"]]
chR <- breaks[["right"]][["54830"]]
lcol <- col.breaks[["left"]][["54830"]]
rcol <- col.breaks[["right"]][["54830"]]


plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["56"]]],
              breaks = breaks, col.breaks = col.breaks, gen="56",
              lg.col = "dodgerblue1", sm.col = "firebrick1",
              region = NULL,
              chrLengths = "./dataFiles/chrLengths_mm10.txt",
              genoFile = "./dataFiles/chr19.filtered.dosage",
              legend = "./dataFiles/chr19.txt.legend.txt",
              outfilePrefix = "./figures/familyHaplos_gen",
              pdfInches = 4)

