### useful for making abe's "chromosome plots" without the chromosome
plot(1:100, rep(15,100), pch='|', col=rainbow(100), ylim=c(10,20))
plot(1:100, rep(15,100), pch='|', col=rainbow(100), ylim=c(10,20), cex=2)
plot(1:100, rep(15,100), pch='|', col=rainbow(100), ylim=c(10,20), cex=5)
points(1:100, rep(13,100), pch='|', col=rev(rainbow(100)), cex=5)

#files <- paste0("/group/palmer-lab/AIL/GBS/dosage/onlyEmpirical/chr", 1:19,
#                ".filtered.snpinfo")
#pos <- read.table("./dataFiles/chr19.filtered.snpinfo", as.is=T, header=T)[3]

# 1. change names of haplo lists to mouse IDs
# 2. extract list of relatives from each element in kincheck.allgenoScaled$idx2
# 3. find out how many relatives there are; make vector 1:length(rel)-1 to plot
# a series of bars for each mouse at regularly spaced intervals.


# read in snp positions, kinCheck
pos <- as.integer(sub("ail.chr19.", "", alt.allele$snp))

# for mouseName in listofRelatedMice
mouse1 <- cbind.data.frame(as.factor(testRun[[1]]), pos)
names(mouse1) <- c("hap", "pos")
split.data <- split(mouse1, f=mouse1$hap)

# make yaxis labels = mouse IDs
plot(split.data$LG$pos/1e6, rep(1, length(split.data$LG$pos)), pch="|",
     xlab = "position (Mb)", ylab="mouse ID", col="blue", cex=2)

# for i in length(relatives)-1
points(split.data$SS$pos/1e6, rep(1, length(split.data$SS$pos)), pch="|",
       col="goldenrod2", cex=2)
points(split.data$LS$pos/1e6, rep(1, length(split.data$LS$pos)), pch="|",
       col="chartreuse3", cex=2)
