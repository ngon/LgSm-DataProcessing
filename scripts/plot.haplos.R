#### FUNCTION: PLOT.HAPLOS -----------------------------------------------------
# plot.haplos takes a list object created by hap.to.color.R and plots haplotypes
# for family members listed in 'ids'.
# You must have ids object loaded to run. which can be an element of sibPairs or
# fams (both are lists of vector objects) which are located in ./pedigree/
# genotypedRelatives.Rdata. Contents of each list are described further in the
# README file in pedigree folder. Alternatively, you could feed it a list of
# any set of IDs. But since the function operates on each vector of ids at a
# time, all elements of the vector are plotted in the same figure.
    # chr is the chromosome name
    # ids (See above)
    # breaks and col.breaks are elements of an object (i named mine "hapdata" as
        # an example) returned by hap.to.color.R.
    # gen is a CHARACTER, e.g. "56". this is if i want to make separate files
        # for each generation. see example at the end.
    # lg.col and sm.col must match the colors used in hap.to.color.R
    # pdfInches = dimensions in inches - the figure is a square.
    # the rest are file paths - self-explanatory. default is chr 19 because it's
        # the smallest chr.

plot.haplos <- function ( chr = 19, ids = sibPairs[[1]], breaks = hapdata[["breaks"]],
                          col.breaks = hapdata[["col.breaks"]], gen=NULL,
                               lg.col="#33a02c", sm.col="#a6cee3",
                               chrLengths = "./dataFiles/chrLengths_mm10.txt",
                               genoFile = "./dataFiles/chr19.filtered.dosage",
                               legend = "./dataFiles/chr19.txt.legend.txt",
                               outfilePrefix = "./figures/familyHaplos_",
                               pdfInches = 4) {

    #### PREPARE DATA #### -----------------------------------------------------
    # get chromosome length
    chrlength <- read.table(chrLengths,as.is=T)$V2[chr] / 1e6

    # get Mb positions of genotyped SNPs in the haplo data set
    # pos$V1 is SNP name to match to typed.snps$V1. pos$V2 = bp position.
    pos <- read.table(legend, header=F)[1:2]
    typed.snps <- read.table(genoFile, header=F, as.is=T)[1]
    pos <- pos[pos$V1 %in% typed.snps$V1,]
    pos <- pos$V2/1e6
    rm(typed.snps) # cleanup

    # Make data frame of all breaks and try to minimize swapping colors from
    # one chromosome to the other.
    all.breaks <- c()
    for (mouse in seq_along(breaks[["left"]])) {
        all.breaks[[mouse]] <- sort(union(breaks[["left"]][[mouse]], breaks[["right"]][[mouse]]))
        all.breaks[[mouse]] <- data.frame(breaks=all.breaks[[mouse]],
                                          left=rep(NA, length(all.breaks[[mouse]])),
                                          right=rep(NA, length(all.breaks[[mouse]])))
        all.breaks[[mouse]][["left"]][match(breaks[["left"]][[mouse]],
                                    all.breaks[[mouse]][["breaks"]])] = col.breaks[["left"]][[mouse]]
        all.breaks[[mouse]][["right"]][match(breaks[["right"]][[mouse]],
                                    all.breaks[[mouse]][["breaks"]])] = col.breaks[["right"]][[mouse]]
    } # for (mouse in seq_along(breaks[["left"]]))
    names(all.breaks) <- names(breaks[["left"]])

    # fill in NA values
    for (mouse in seq_along(all.breaks)) {
        for (i in seq_along(all.breaks[[mouse]][["left"]])){
            if(is.na(all.breaks[[mouse]][["left"]][i])) {
                all.breaks[[mouse]][["left"]][i] <- all.breaks[[mouse]][["left"]][i-1]
            } # if(is.na(all.breaks[[mouse]][["left"]])
        } # for (i in seq_along(all.breaks[[mouse]][["left"]]))
        for (i in seq_along(all.breaks[[mouse]][["right"]])){
            if(is.na(all.breaks[[mouse]][["right"]][i])) {
                all.breaks[[mouse]][["right"]][i] <- all.breaks[[mouse]][["right"]][i-1]
            }# if(is.na(all.breaks[[mouse]][["right"]])
        }# for (i in seq_along(all.breaks[[mouse]][["right"]]))
    } # for (mouse in seq_along(all.breaks))


    #### PLOT HAPLOS #### ------------------------------------------------------
    # initialize PDF
    pdf(file=paste0(outfilePrefix, gen, ".pdf"), width=pdfInches, height=pdfInches)

    for (fam in seq_along(ids)) {
        # extract info from all.breaks one family at a time
        sibtest <- c(all.breaks[names(all.breaks) %in% ids[[fam]]])

        # set up plot area.
        par(font = 2, font.axis = 2, font.lab = 2, las = 1,
            plt = c(0.2, 0.9, 0.12, 0.95), mar=c(2,4,1,2.7)+0.1)

        # draw chromosome skeletons and label the y-axis
        plot(1, 1, col = 0, xlim = c(1, (2 * length(sibtest) + 1)),
             ylim = c(0, chrlength+2), xaxt = "n", yaxs = "i", xlab = "", ### VAR chrlength
             ylab = paste0("Chr ", chr, " position (Mb)"), cex=0.9) ### VAR chr
        for (a in 1:length(sibtest)) {
            lines(2 * c(a, a), c(0, chrlength))
        } # for (a in 1:length(sibtest))

        # gridlines
        abline(h = 0:20 * 5, col = "grey80")

        # add a 'legend' (if you use the legend function in baseR, it takes up
        # a lot of space and the position may need to be adjusted for different
        # sizes of families. putting labels in the margin avoids this issue).
        mtext(side = 4, at = c((chrlength/2)+5, chrlength/2), adj = -0.2,
              text = c("LG/J", "SM/J"), col = c(lg.col, sm.col), cex = 0.8)

        # x-axis labels for each mouse in sibtest
        for (mouse in seq_along(names(sibtest))) {
            mtext(side = 1, at = mouse * 2, line = 0.2, text = names(sibtest)[mouse])
        } # for (mouse in seq_along(names(sibtest)))

        for (mouse in seq_along(sibtest)) {
            offset = 0.5
            xm = 2 * mouse
            xl = xm - offset
            xr = xm + offset

            # plot the first block on each side of the chromosome
            rect(xl, pos[sibtest[[mouse]][["breaks"]][1] ], xm,
                 pos[sibtest[[mouse]][["breaks"]][2] ],
                 col=sibtest[[mouse]][["left"]][1],
                 density=NA, border=NA)
            rect(xm, pos[sibtest[[mouse]][["breaks"]][1] ], xr,
                 pos[sibtest[[mouse]][["breaks"]][2] ],
                 col=sibtest[[mouse]][["right"]][1],
                 density=NA, border=NA)

            # go through each break and see if the next color should be swapped
            for (i in 2:nrow(sibtest[[mouse]])) {
                if (sibtest[[mouse]][["left"]][i] != sibtest[[mouse]][["right"]][i]) {
                    if(sibtest[[mouse]][["left"]][i] == sibtest[[mouse]][["right"]][i-1] |
                       sibtest[[mouse]][["left"]][i-1] == sibtest[[mouse]][["right"]][i]){
                        tmp <- sibtest[[mouse]][["left"]][i]
                        sibtest[[mouse]][["left"]][i] <- sibtest[[mouse]][["right"]][i]
                        sibtest[[mouse]][["right"]][i] <- tmp
                    } #  if sibtest ..... | .....
                } # if (sibtest[[mouse]][["left"]][i] != .....
                rect(xl, pos[sibtest[[mouse]][["breaks"]][i] ], xm,
                     pos[sibtest[[mouse]][["breaks"]][i+1] ],
                     col=sibtest[[mouse]][["left"]][i],
                     density=NA, border=NA)
                rect(xm, pos[sibtest[[mouse]][["breaks"]][i] ], xr,
                     pos[sibtest[[mouse]][["breaks"]][i+1] ],
                     col=sibtest[[mouse]][["right"]][i],
                     density=NA, border=NA)
            } # for (i in 2:length(sibtest[[mouse]][["breaks"]]))
        } # for (mouse in seq_along(sibtest))

    } # for (fam in seq_along(ids))
    dev.off()
} # plot.haplos - END OF FUNCTION

###################### EXAMPLES FOR RUNNING PLOT.HAPLOS ########################

## RUN HAP.TO.COLOR
# source("./scripts/hap.to.color.R")
# load("./dataFiles/chr19.hap.RData")
# hapdata <- hap.to.color(hap=chr19.hap)
# names(hapdata) <- c("breaks", "col.breaks")
# rm(chr19.hap)

############## GENERATIONS 50 - 56 (SIB PAIR DATA) ############################
# i want to put the plots into several different PDF files to avoid making them
# too large to open. here i'm grouping them by generation. here i borrowed code
# from ./pedigree/genotyped.relatives.R: i load the pedigree and get info for
# mice that were genotyped and have siblings (e.g. mice in sibPairs). I make
# sure that the list of dams in info is identical to names(sibPairs). Now I can
# determine what range of sibPair elements comprises each generation of AILs
# and plot them in separate PDF files.
#
# load("./pedigree/info.phenoSamples.RData")
# load("./pedigree/genotypedRelatives.RData")
#
# info <- info[order(info$gen, info$dam),]
# ids <- read.table("./genotyped.samples.txt", as.is=T)$V1
# info$id <- info$id - 0.1
# info <- info[info$id %in% ids,]
# test <- info[info$dam %in% as.numeric(names(sibPairs)),c(1:3)]
# test <- test[which(!duplicated(test$dam)),]
# test <- test[order(test$dam),]
# identical(test$dam, as.integer(names(sibPairs))) # yes
# test$gen <- as.factor(test$gen)
# row.names(test) <- 1:509
# test <- split(test, test$gen, drop=T)
# genrows <- lapply(test, row.names)
# genrows <- lapply(genrows, as.integer)
# rm(test, info, ids)

# plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["56"]]],
#               gen = "56",pdfInches = 5)
# plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["56"]]],
#               gen = "56",pdfInches = 5)
# plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["55"]]],
#               gen = "55",pdfInches = 5)
# plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["54"]]],
#               gen = "54",pdfInches = 5)
# plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["53"]]],
#               gen = "53",pdfInches = 5)
# plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["52"]]],
#               gen = "52",pdfInches = 5)
# plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["51"]]],
#               gen = "51",pdfInches = 5)
# plot.haplos ( chr = 19, ids = sibPairs[seq_along(sibPairs) %in% genrows[["50"]]],
#               gen = "50",pdfInches = 5)

############## GENERATIONS 39 - 43 (SIB PAIR DATA) ############################
# This time I'll group the families together by the parents' generation. There
# will be a separate PDF for parents born in F39 - 42.

# ped<- read.table("./pedigree/pedforQTLRel.txt", header=T)

# the dam is the last entry in each fams vector. get her id, then her gen.
# mamaMouse <- c()
# for (family in seq_along(fams)){
#     mamaMouse <- c(mamaMouse, fams[[family]][length(fams[[family]])])
# }
# mamaMouse <- mamaMouse + 0.1
# mmGen <- ped[ped$id %in% mamaMouse,]
# mmGen <- mmGen[c(1,5)]
# identical(mmGen$id, mamaMouse)
# test <- sort(mamaMouse, decreasing=FALSE)
# identical(test, mamaMouse) # yes, mamaMouse is in decreasing order
# mmGen <- mmGen[order(mmGen$id),]
# identical(mmGen$id, mamaMouse) # now TRUE
# mmGen$id <- mmGen$id - 0.1
# names(fams) <- mmGen$id
    # SAVING FAMS WITH DAMS AS LIST NAMES FOR FUTURE USE
    # save(fams, sibPairs, file="./pedigree/genotypedRelatives.RData")
# row.names(mmGen) <- 1:190
# mmGenSplit <- split(mmGen, mmGen$generation, drop=T)
# genrows <- lapply(mmGenSplit, row.names)
# genrows <- lapply(genrows, as.integer)
# rm(mamaMouse, mmGen, test, mmGenSplit, ped)

# plot.haplos ( chr = 19, ids = fams[seq_along(fams) %in% genrows[["F39"]]],
#               gen = "39", pdfInches = 7)
# plot.haplos ( chr = 19, ids = fams[seq_along(fams) %in% genrows[["F40"]]],
#               gen = "40", pdfInches = 6)
# plot.haplos ( chr = 19, ids = fams[seq_along(fams) %in% genrows[["F41"]]],
#               gen = "41", pdfInches = 6)
# plot.haplos ( chr = 19, ids = fams[seq_along(fams) %in% genrows[["F42"]]],
#               gen = "42", pdfInches = 6)
