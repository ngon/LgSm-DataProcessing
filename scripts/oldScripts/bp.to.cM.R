### Purpose: convert bp to cM distance and save info
positions <- c()
chromosomes <- paste0("chr", 1:19)
fileNames <- c()

# prepare to read each chr's file
for (i in seq_along(chromosomes)){
    fileNames[i] <- paste0(chromosomes[i], ".filtered.snpinfo")
    positions[i] <- read.table(fileNames[i], as.is=T, header=F, sep="\t")[2]
}

# give each element of the list a name and convert to df with an additional
# column for chromosome number.
names(positions)[1:19] <- chromosomes
positions <- lapply(positions, data.frame)

for (i in seq_along(names(positions))){
    positions[[i]][[2]] <- rep(i, times=length(positions[[i]][[1]]))
}

# switch the order of columns so that chr is col 1 and pos is col 2
positions <- lapply(positions, function(data) { data <- data[c(2,1)]} )

# write to 19 separate files
for (i in seq_along(positions)){
write.table(positions[[i]], file=paste0("chr", i, "bp.txt"), sep="\t",
            col.names=F, row.names=F, quote=F)
}

# NEXT i simply loaded each file onto the JAX mouse map converter at
# http://cgd.jax.org/mousemapconverter/. I converted bp position (GRCm38)
# to Sex-averaged cM (Cox, 2014). I downloaded the conversion file for each
# chromosome and merged them using the commands below.
# New files are located in CRI within /AIL/GBS/dosage/onlyEmpirical/centimorgams

snpinfo <- c()
for ( chr in 1:19 ){
bp <- read.table(file=paste0("./chr", chr, "bp.txt"), sep="\t", header=F, as.is=T)
cm <- read.table(file=paste0("./chr", chr, ".cox.txt"), sep="\t", header=F, as.is=T)[2]
snpinfo[[chr]] <- data.frame(bp, cm)
}
names(snpinfo)[1:19] <- chromosomes

for (i in seq_along(snpinfo)){
    write.table(snpinfo[[i]], file=paste0("chr", i, ".bp.cM.txt"), sep="\t",
                col.names=F, row.names=F, quote=F)
}

for (i in seq_along(snpinfo)){
    names(snpinfo[[i]])[1:3] <- c("chr", "bp", "cM")
}

save(snpinfo, file="snpChrPosCm.RData")





