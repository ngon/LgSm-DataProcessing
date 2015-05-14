### Purpose: convert bp to cM distance and save info

chromosomes <- paste0("chr", 1:19)
fileNames <- c()
for (i in seq_along(chromosomes)){
    fileNames[i] <- paste0(chromosomes[i], ".filtered.snpinfo")
}


get.snpinfo <- function(filename) {

    positions[[chr]] <- read.table(filename, as.is=TRUE, sep="/t")[2]




}

