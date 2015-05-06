cmds <- c()

choose.traits <- function(traits=list(), filename){
   chosen.covs <- c()
   index.pheno <- c()

   for (trait in traits){
        covs <- unlist(traitcovs[[trait]])
        chosen.covs <- append(chosen.covs, covs)
        index.pheno <- append(index.pheno, which(pheno.names == trait))
        }
   chosen.covs2 <- chosen.covs[-(which(duplicated(chosen.covs)))]

   cov.df <- covars[, names(covars)%in%chosen.covs2]

   write.table(cov.df, file=paste0(filename, ".emp.covs"),
                sep="\t", quote=F, row.names=F, col.names=F)

   for (chrom in 1:19) {
       cmds <- c(cmds, paste0("gemma -g /group/palmer-lab/AIL/GBS/dosage/chr", chrom,
                              ".filtered.dosage.emp -p /group/palmer-lab/AIL/LgSm-DataProcessing/phenos.allgeno.txt -k /group/palmer-lab/AIL/qtlmapping/kinship/chrNot",
                              chrom,".cXX.txt -a /group/palmer-lab/AIL/GBS/dosage/chr",
                              chrom, ".filtered.snpinfo.emp -c /group/palmer-lab/AIL/qtlmapping/covariates/",
                              trait, ".emp.covs -lmm 2 -maf ", MAF, " -o ", trait, ".chr",
                              chrom, " -n ", index.pheno[1:length(index.pheno)]))
   }

write.table(cmds, file=paste0("./gemma.multivar.emp.cmds"),
            row.names=F, col.names=F, quote=F)

}


choose.traits(traits=c("tail", "habituation", "glucose"), filename="tail_hab_glu")

