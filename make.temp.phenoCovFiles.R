## May 1 2015 ##

## Summary ##
## I'm updating phenotype and covariate files in LgSm-DataProcessing so that I
## can quickly map QTL on only empirical SNPs over the weekend. I'm not
## documenting everything I've done in detail yet because I need to do some
## regression tests and make residual plots for covariates for all of the --.age
## covariates. Why? Because there are some older animals and younger animals
## that should be excluded from the phenotype data frame. Older mice will need
## to be thrown out for CPP for sure, and probably glucose/weight I'm less sure
## if age will have a large effect on PPI, but I want to test all traits for
## an age effect empirically. I don't have time to do it right away, so I'm
## excluding mice at the extreme ends of the age distribution.

## NOTE: I am making changes IN THE TOP LEVEL DIRECTORY ONLY, not in dataFiles.
## Final changes will be made after I get GEMMA running, accompanied with a
## detailed readme file explaining changes.

pheno <- read.table("phenotypes.orm.txt", header=T, as.is=T)
cov <- read.table("./dataFiles/phenotypes/covariates.txt", header=T, as.is=T)
quantile(cov$cpp.age, probs=seq(0,1,0.05), na.rm=T)

covs <- cov[cov$cpp.age <47,]
covs2 <- covs[covs$cpp.age <=91,]
ids <- c(covs$id, covs2$id)
ids <- sort(ids)

cov.new <- cov[!cov$id %in% ids,]
pheno.new <- pheno[!pheno$id %in% ids,]

write.table(cov.new, file="covariates.txt", quote=F, sep="\t", row.names=F,
            col.names=TRUE)
write.table(pheno.new, file="phenotypes.orm.txt", quote=F, sep="\t", row.names=F,
            col.names=TRUE)



