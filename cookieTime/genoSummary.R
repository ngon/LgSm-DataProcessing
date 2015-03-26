# IMPUTED VS EMPIRICAL SNPS --------------------------------------
# MAF and HZ for each
# /AIL/GBS/preimpute/file.geno  <- empirical
# /AIL/GBS/dodage/file.filtered.dosage <- imputed and empirical

#  MAF ------------------------------------------------------------
# Returns the minor allele frequency given a vector of genotypes
# encoded as allele counts.
computemaf <- function (geno) {
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f,1-f))
}
# HETEROZYGOSITY BY SAMPLE ----------------------------------------
# use a threshold of 0.8-1.2 for hets

# HETEROZYGOSITY BY SNP -------------------------------------------
# use a threshold of 0.8-1.2 for hets

# PED VS GRM RELATEDNESS ------------------------------------------
# matrix(snp, sample) : cor(r) between cols(samples) = GRM correl matrix

# LD DECAY -------------------------------------------------------
# transpose GRM correl matrix to get cor(r) between snps 
# AIL/knownSNPs/imputeHaplotypes
# .hap files, 0=ref, 1=nonref
# .legend files give row names (snps)






