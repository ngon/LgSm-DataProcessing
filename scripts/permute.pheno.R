### PERMUTE PHENOTYPES ###
### PURPOSE: Create files with permuted phenotypes and covariates to feed
### GEMMA. Goal is to find significance thresholds for GWAS data to include
### in my 2015 CTC poster. Five GWAS will be run for each phenotype.I will
### derive a threshold based on the pvalues that fall in the 95th percentile
### in permuted data sets.

### LOAD DATA ------------------------------------------------------------------
pheno <- read.table("./phenotypes.orm.txt", sep="\t", header=TRUE, as.is=TRUE)
cov <- read.table("./covariates.orm.txt", sep="\t", header=TRUE, as.is=TRUE)
load("./traitcovs.RData")

traits.to.perm <- c("act1.t", "act2.t", "act4.t", "act5.t", "act8.t", "cpp.diff",
                    "sc8.t", "sc1.t", "startle", "wild.binary", "ppi6.logit",
                    "ppi12.logit", "tail","glucose")

# pheno data for traits.to.perm
phenew <- which(names(pheno) %in% names(traits.to.perm))
phenew <- data.frame(pheno$id, pheno[phenew])

# one list of covariates for each trait in traits.to.perm
tnames <- traitcovs[names(traitcovs) %in% traits.to.perm]

# list of dfs, each with a col for the trait and its covariates
phco <- list()
for (trait in traits.to.perm){
    get.covs <- cov[,which(names(cov) %in% tnames[[trait]])]
    phco[[trait]] <- data.frame(phenew[[trait]], get.covs)
    names(phco[[trait]])[1] <- trait
}

# for each column, define 'NA' and data entries as separate blocks
# plots will be each col in a df; blocks will be NA and data.
# permute within blocks and within plots.

# outputs row names of NA and DATA entries in separate vectors
get.blocks <- function(trait){
    blks <- c()
    blockNA <- c()
    blockData <- c()
    for (col in seq_along(phco[[trait]])){
        blockNA[[col]] <- which(is.na(phco[[trait]][col]))
        blockData[[col]] <- which(!is.na(phco[[trait]][col]))
        blks[[trait]] <- c(blks[[trait]],list(blockNA[[col]], blockData[[col]]))
    }
    return(blks)
}
# within each trait, odd blocks are NA and even blocks are data
# access as  > blocks[['startle']][3]
blocks <- sapply(traits.to.perm, get.blocks)
names(blocks)[1:14] <- traits.to.perm

# outputs a numeric logical vector for each column
get.levels <- function(trait){
    lvs <- c()
    for(col in seq_along(phco[[trait]])){
        lvs[[col]] <- as.integer(is.na(phco[[trait]][col]))
    }
    return(lvs)
}
# access same as blocks, e.g. > levs[['startle']][1]
levs <- sapply(traits.to.perm, get.levels)

# permutations with permute package
library(permute)

# test case on 1 column of data: it works.
# foo.data <- phco[["startle"]]
# foo.rows <- blocks[["startle"]]
#
# blockS <- c(levs[['startle']][1], levs[['startle']][2])
#
# ctrl <- how(within=Within(type="free", constant=TRUE),
#             blocks=as.factor(blockS[[1]]), nperm=5)
#
# a<- shuffle(n=1096, control=ctrl)
# b<- shuffle(n=1096, control=ctrl)
#
# test <- data.frame(blockS[[1]], a, b)
# names(test)[1] <- "c"
# check<- test[test$c == 1,]
# check <- check[do.call(order, check),]
# aa<- sort(check$a)
# bb<- sort(check$b)
# cbind(aa,bb, blocks[[startle]][1])


# FUNCTION: PERMUTE.PHENOS ----------------------------------------------------
# permute phenotypes, leaving NA entries constant.

permute.phenos <- function(trait, data, levels, nset=5){
    require(permute)
    thetrait <- phco[[trait]][[1]]
    thecovars <- phco[[trait]][-1]
    thelevels <- levs[[trait]][[1]]
    n <- length(thelevels)

    # permute rows, keeping NA values constant. the output, perms, is a list
    # of nset permutations for each column of trait data
    ctrl <- c()
    perms <- c()
        ctrl[[trait]] <- how(within=Within(type="free", constant=TRUE),
                    blocks=as.factor(thelevels))
        perms[[trait]] <- shuffleSet(n, nset=nset, ctrl[[trait]])

    p.perm <- c()
    c.perm <- c()
    for(i in 1:nset){
        p.perm[[trait]][[i]] <- phco[[trait]][[trait]][(perms[[trait]][i,])]
        c.perm[[trait]][[i]] <- thecovars[(perms[[trait]][i,]),]
    }
    p.perm[[trait]] <- as.data.frame(p.perm[[trait]])
    list(p.perm, c.perm)
    }




### PERMUTE PHENOS AND SAVE DATA ----------------------------------------------
permData <- sapply(traits.to.perm, permute.phenos, data=phco,
                       levels=levs, nset=5 )

# separate phenotypes from covariates
odds <- which(seq_along(permData) %% 2 !=0)
evens <- which(seq_along(permData) %% 2 == 0)

covp <- permData[c(evens)]
phenop <- permData[c(odds)]

# change pheno list to data frame, rename columns, and save as a text file.
phenop <- as.data.frame(phenop)
nms<-c()
for (i in traits.to.perm){ nms <- c(nms, rep(i, times=5)) }
pnames <- paste0("perm", 1:5, ".", nms)
names(phenop) <-pnames
write.table(phenop, file="./permutedPhenotypes.txt", sep="\t",
            row.names=F, col.names=T, quote=F)

# make a separate file for each set of permuted covariates
for(j in seq_along(1:5)){
for(i in seq_along(covp)){
    cov <- covp[[i]][[1]][[j]]
    write.table(cov, file=paste0("./perm", j, ".", traits.to.perm[i], ".cov.txt"),
                row.names=F, col.names=T, quote=F)
}
}

### base R perms - unused -----------------------------------------------------
# permute.pheno <- function(phenotype, phenoData, covData){
#
#     phenoData <- phenew
#     data <- data.frame(phenew[[phenotype]])
#     names(data)[1] <- phenotype
#     covData <- traitcov[[phenotype]]
#
#     a<- rownames(phenew[which(!phenew[['wild.binary']] == "NA"),])
#
#     sample(phenew[['wild.binary']][which(!phenoData[['wild.binary']] == "NA")])
#
#     samp <- list()
#     perms <- list()
#     p.names <- c()
#
#     for (i in 1:5){
#         p.names[i] <- paste0(i, phenotype)
#         p.samp[[i]] <- sample(phenoData[[phenotype]][which(!phenoData[[phenotype]] == "NA")])
#         p.perm[[i]] <- replace(x=phenoData[[phenotype]],
#                               list=which(!phenoData[[phenotype]] == "NA"),
#                               values=samp[[i]])
#
#         c.samp[[i]] <- order()
#         c.perm[[i]] <- replace()
#
#     }
#
#     covData <- traitcov[[phenotype]]
#
#     names(p.perm) <- p.names
#     list(p.perms, c.perm)
# }
#
# testPerm <-lapply(traits, permute.pheno)


### GET P VALUES FROM GWAS ON PERMUTED PHENOS ----------------------------------

### FUNCTION: get.filenames ---------------------------------------------------
### (get permed assoc files to feed to getP)
get.filenames <- function(trait){
    filenames <- c()
    for (i in 1:5){
        filenames[i] <- paste0("/group/palmer-lab/AIL/qtlmapping/output/permutedTraits/perm",
                               i, ".")
    }
    return(filenames)
}
files <- lapply(traits.to.perm, get.filenames)
files <- unlist(files)


### FUNCTION: get.P ------------------------------------------------------------
### (adapted from getP.plotQ.R - this one does not plot qqs)

getP.plotQ <- function(directory=getwd(), trait){
    # match file names to trait
    fileStart <- c()
    fileNames <- c()
    chromosomes <- 1:19

    for (j in seq_along(1:5)){
        for(i in seq_along(chromosomes)){
            fileStart[j] <- paste0("/group/palmer-lab/AIL/qtlmapping/output/permutedTraits/perm",
                                   j, ".")
            fileNames[i] <- paste0(fileStart[j], trait, ".chr", chromosomes[i],
                                   ".assoc.txt")
        }
    }
    # read only p_lrt columns from selected files and compile them into a single data frame
    chosenFiles <- lapply(file.path(directory, fileNames), read.table, sep="\t",
                          header=T, colClasses = c(rep("NULL", 8), "numeric"))
    pValues <- do.call(what=rbind.data.frame, args=chosenFiles)

    # take the -log10 and get expected and observed values
    pValues <- pValues$p_lrt
    obs = -log10(sort(pValues,decreasing=F))
    exp = -log10( 1:length(obs)/length(obs) )

}














