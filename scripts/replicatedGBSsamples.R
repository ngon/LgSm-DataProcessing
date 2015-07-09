ped <- read.table("./pedigree/pedforQTLRel.txt", header=T)
ids <- read.table("./dataFiles/replicateGenos.txt", header=F)[1]
ids<- as.character(ids)

test <- (lapply(ids, strsplit, split=".", fixed=T))

samples <- c()
for (i in seq_along(test)){
    samples <- c(samples, test[[i]][[1]][1])
}

id.tmp <- unique(samples)

# save control animals and regenotyped mice in separate structures
controls <- data.frame(id = c(45905, 45915, 45943, 45946, 26874, 46679),
                       gen = c( 'LG', 'LG', 'SM', 'SM', 'F1', 'F1'))

id.tmp[which(nchar(id.tmp) > 5)] <- rep(0, 9)
id.tmp <- sort(as.numeric(id.tmp))
id.tmp <- id.tmp[10:length(id.tmp)]
ids <- id.tmp + 0.1
rm(id.tmp, test)

rep.info <- ped[ped$id %in% ids,]
rep.info$id <- rep.info$id - 0.1
rep.info <- merge(rep.info, controls, all=TRUE)

