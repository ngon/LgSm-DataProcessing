### match pheno and genoIDs


# read data
pheno <- read.table("./ailPhenotypes.txt", sep="\t", header=FALSE, na.strings="NA")
phenoNames <- read.table("./phenColNames.txt", sep="\t")
names(pheno) <- phenoNames[1:100,1]

geno <- read.table("./genoIDs.txt", sep="", header=FALSE)
geno <- t(geno)
geno <- as.data.frame(geno)
names(geno)[1] <- "id"

# match ids in pheno and genoID files
phenoMatchGeno <- merge(geno, pheno, by="id", all=TRUE)
phenoMatchGeno <- phenoMatchGeno[,2:100]

# create new files 
write.table(phenoMatchGeno, "./phenoMatchGeno.txt", sep="\t",
            row.names=FALSE, quote=FALSE)

write.table(phenoMatchGeno$id, "./idsPhenoMatchedGeno.txt", sep="\t")

