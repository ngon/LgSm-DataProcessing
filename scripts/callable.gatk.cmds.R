# first run list.files(pattern=pattern) to check filenames for the mouse with reps
# count the number of reps and change values in the variable 'one', in names of
# df.dups columns, and in the cmds for loop. .

one <- rep(paste0("rep", 0:15), times=16)
tmp <- data.frame(one, sort(one))
names(tmp)[2] <- "two"

dups <- c()
for (j in unique(tmp$two)) {
    dups[[j]] <- which(tmp$one == j)
}
df.dups <- as.data.frame(dups)

df.dups <- df.dups[c("rep0", "rep1", "rep2", "rep3", "rep4", "rep5", "rep6", "rep7",
                     "rep8", "rep9", "rep10", "rep11", "rep12", "rep13", "rep14","rep15")]

#df.dups<- (df.dups[c(2,13,18:24,3:12,14:17,1)])
#df.dups <- df.dups[c(24,1:23)]
m.dups<- as.matrix(df.dups)
m.dups[upper.tri(m.dups, diag=TRUE)] <- NA

pairs <- tmp[(which(!is.na(m.dups))),]

cmds <- c()

for (i in 1:length(pairs$one)){
    cmds <- append(cmds, paste0("java -Xmx4g -jar /apps/software/GenomeAnalysisTK/3.3-0/GenomeAnalysisTK.jar -R /group/palmer-lab/reference_genomes/mouse/mm10.fasta -T CompareCallableLoci -comp1 /group/palmer-lab/AIL/GBS/callableLoci/26874.", pairs[i,1],".isCallable.bed -comp2 /group/palmer-lab/AIL/GBS/callableLoci/26874.", pairs[i,2], ".isCallable.bed -L chr19 -o 26874_", pairs[i,1], "_", pairs[i,2], ".table"))

}

