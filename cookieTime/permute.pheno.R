pheno <- read.table("./phenotypes.orm.txt", sep="\t", header=TRUE, as.is=TRUE)
cov <- read.table("./covariates.orm.txt", sep="\t", header=TRUE, as.is=TRUE)
load("./traitcovs.RData")

traits.to.perm <- c("act1.t", "act2.t", "act4.t", "act5.t", "act8.t", "cpp.diff", "sc8.t",
            "sc1.t", "startle", "wild.binary", "ppi6.logit", "ppi12.logit", "tail",
            "glucose")

tnames <- traitcovs[names(traitcovs) %in% traits.to.perm]

phenew <- which(names(pheno) %in% names(traits.to.perm))
phenew <- data.frame(pheno$id, pheno[phenew])



# nalist <- list()
# for (name in names(phenew)){
#     nalist[[name]] <- which(is.na(phenew[[name]]))
# }


permute.pheno <- function(phenotype, phenoData, covData){

    phenoData <- phenew
    covData <- traitcov[[phenotype]]
    #data <- read.table(file, sep="\t", header=T, as.is=T)

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







