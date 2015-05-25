pheno <- read.table("./phenotypes.orm.txt", sep="\t", header=TRUE, as.is=TRUE)
cov <- read.table("./covariates.orm.txt", sep="\t", header=TRUE, as.is=TRUE)
load("./traitcovs.RData")

traits.to.perm <- c("act1.t", "act2.t", "act4.t", "act5.t", "act8.t", "cpp.diff",
                    "sc8.t", "sc1.t", "startle", "wild.binary", "ppi6.logit",
                    "ppi12.logit", "tail","glucose")

# pheno data for traits.to.perm
phenew <- which(names(pheno) %in% names(traits.to.perm))
phenew <- data.frame(pheno$id, pheno[phenew])

# nalist <- list()
# for (name in names(phenew)){
#     nalist[[name]] <- which(is.na(phenew[[name]]))
# }

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

[[i]][[1]][i]
# make a separate file for each set of permuted covariates

for(i in seq_along(covp)){
    cov <- covp[[i]][[1]][[5]]
    write.table(cov, file=paste0("./perm5.",traits.to.perm[i], ".cov.txt"),
                row.names=F, col.names=T, quote=F)
}


### PERMUTE COVARS AND SAVE DATA ----------------------------------------------
permedCovars <- sapply(traits.to.perm, permute.phenos, data=phco,
                       levels=levs, nset=5 )
names(permedCovars) <- traits.to.perm
covarSet <- permedCovars[[trait]][[1]]


# returns a list with nset rows and n columns for each trait vector.
# level 1: for each phenotype in traits.to.perm
# level 2: for each column of data assoc w/ traits in traits.to.perm
# level 3: matrix of row indexes for each cov or pheno column in phco
# > test[[1]][1,]
#test <- permute.phenos(trait="startle", data=phco, levels=levs, nset=5)
permutedRows <- sapply(traits.to.perm, permute.phenos, data=phco,
                       levels=levs, nset=5)
# name structures for easier indexing
for (trait in traits.to.perm){
    names(permutedRows[[trait]]) <- names(phco[[trait]])
}
# access a set of permuted data for a given trait/cov as so:
permutedRows[['act1.t']][['act1.t']][1,]


get.permedData <- function(trait, nset){
p.perm <- c()
for(i in 1:nset){ # TO DO: replace permutedRows with perms
    p.perm[[trait]][[i]] <- phco[[trait]][[trait]][(permutedRows[[trait]][[trait]][i,])]
}
p.perm[[trait]] <- as.data.frame(p.perm[[trait]])
names(p.perm[[trait]])[1:nset] <- paste0("p",trait, 1:nset)
return(p.perm)
}


c.perm <- c()
for(j in names(phco[[trait]][-1])){
    c.perm[[j]] <- cbind(c, phco[[trait]][[j]][(permutedRows[[trait]][[j]][i,])])
}



### base R perms - in progress
permute.pheno <- function(phenotype, phenoData, covData){

    phenoData <- phenew
    data <- data.frame(phenew[[phenotype]])
    names(data)[1] <- phenotype
    covData <- traitcov[[phenotype]]

    a<- rownames(phenew[which(!phenew[['wild.binary']] == "NA"),])

    sample(phenew[['wild.binary']][which(!phenoData[['wild.binary']] == "NA")])

    samp <- list()
    perms <- list()
    p.names <- c()

    for (i in 1:5){
        p.names[i] <- paste0(i, phenotype)
        p.samp[[i]] <- sample(phenoData[[phenotype]][which(!phenoData[[phenotype]] == "NA")])
        p.perm[[i]] <- replace(x=phenoData[[phenotype]],
                              list=which(!phenoData[[phenotype]] == "NA"),
                              values=samp[[i]])

        c.samp[[i]] <- order()
        c.perm[[i]] <- replace()

    }

    covData <- traitcov[[phenotype]]

    names(p.perm) <- p.names
    list(p.perms, c.perm)
}

testPerm <-lapply(traits, permute.pheno)







