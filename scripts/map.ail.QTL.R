# process phenotype for Natalia #

###### CONSTANTS
MAF=0.05

#### Read phenotype file
pheno     <- read.table('phenosAIL.txt')
pheno.names <- read.table("phenColNames.txt", skip=1, as.is=T)[,2]
names(pheno) <- pheno.names

#### Read master covariate file
covars <- read.table("covariatesAIL.txt", sep="\t", as.is=T)
covar.names <- read.table("covColNames.txt", skip=1, as.is=T)[,2]
names(covars) <- covar.names


####### Test code #########
cpp8.t <- pheno$cpp8.t
batch <- as.factor(covars$batch)

anova.cpp <- anova(lm(cpp8.t~batch))

batch13 <- as.factor(covars$batch == 13)
anova.cpp <- anova(lm(cpp8.t~batch13))

#################### Run Gemma for chosen trait ##############
group.name <- "cpp8"
chosen.pheno <- c("cpp8.t")
chosen.covars <- c("one", "sex") #,  paste0("box", 1:11))

index.pheno <- sapply(chosen.pheno, function(x) { which(pheno.names == x)} )
## Make covariate file.
chosen.covars <- sapply(chosen.covars, function(x) { which(covar.names == x)} )
chosen.covars <- covars[, chosen.covars]
write.table(chosen.covars, file=paste0(group.name, ".covs"), sep="\t", quote=F, row.names=F, col.names=F)

## RUN GEemma chromosome wise
cmds <- c("#! /bin/bash\n#$ -cwd\n#$ -j y\n")
for (index in 1:length(index.pheno)) {
    for (chrom in 1:19) {
        cmds <- c(cmds, paste0("/home/shyamg/bin/gemma -g genotypes/ail.chr", chrom, ".dosage -p phenosAIL.txt -k kinship/notChr", chrom,".cXX.txt -a snpinfo/ail.chr", chrom, ".snpinfo -c ", group.name, ".covs -lmm 2 -maf", MAF, " -o ", chosen.pheno[index], ".chr", chrom, ".out -n ", index.pheno[index]))
    }
    write.table(cmds, file=paste0("scripts/", group.name, ".", chosen.pheno[index], ".sh"))
}






#Code to plot manhattanplot
chrLens    <- read.table('/home/shyamg/projects/Mouse/CFW/maps/chrLengths_mm9.txt')$V2
chrLens    <- c(0, cumsum(chrLens*1.0))
labPos     <- (read.table('/home/shyamg/projects/Mouse/CFW/maps/chrLengths_mm9.txt')$V2)/2.0
labPos     <- chrLens[1:19]+labPos

plotManhattan <- function(gemmaout, trait) {
  output <- read.table(gemmaout, header=T)
  pvals <- output$p_lrt
  pvals <- -log10(pvals)
  positions <- chrLens[output$chr]+output$ps
  cols <- output$chr
  odds <- which(cols%%2 == 1)
  cols[which(cols%%2 == 0)] = 'dark red'
  cols[odds] = 'dark blue'
  plot(positions, pvals, col=cols, pch=19, cex=0.5, main=trait, ylab='-log10(p-value)', xlab='Position', axes=F, frame.plot=T)
#  abline(h=-log10(0.05/length(pvals)), lty=2, col='orangered')
  axis(2)
  axis(1, at=labPos, labels=1:19, tick=F)
}

