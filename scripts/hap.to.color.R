
######  FUNCTION: HAP TO COLOR -------------------------------------------------
# A haplotype list must be preloaded before running hap.to.color. I have stored
# these as RData files for each chromosome. Name format is "chrN.hap.RData" and
# the object within is always "chrN.hap". It returns two lists: breaks, with
# list positions where a haplotype switch occurs, and col.breaks, where breaks
# is matched to colors representing LG and SM.
    # hap -- a list of haplotypes (.RData files for each chr are located in
        # /group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/lgsmCrossovers)
    # lg.col -- default is a deep green
    # sm.col -- default is a pale blue

hap.to.color <- function (hap, lg.col="#33a02c", sm.col="#a6cee3") {

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

# Example use:
    # load("./dataFiles/chr19.hap.RData")
    # hapdata <- hap.to.color(hap=chr19.hap)
    # names(hapdata) <- c("breaks", "col.breaks")
    # rm(chr19.hap) # to save memory
# str of the resulting list is: listName [[breaks]] [["right"]] [[mouseID]]
