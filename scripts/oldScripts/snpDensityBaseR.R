## Plot SNP density by chromosome in base R

ChrDensity <- function(ChrLength, Positions, Values, ChromosomeBgColor="gray90", ChromosomeBorderColor="black", ValuesColors="black", PlottingGenes=F, GeneColors="black", ylab="Density", xlab="Chromosomal position (Mb)", ...) {

    PositionsMb <- Positions / 1000000
    ChrLengthMb <- ChrLength / 1000000

    plot(NA,NA, xlab=xlab, ylab=ylab , bty="n", xlim=c(0,ChrLengthMb), ylim=c(-1.1,max(Values)), ...)
    rect(0, -0.1, ChrLengthMb,-1, col=ChromosomeBgColor, border=ChromosomeBorderColor)
    nbG <- length(Positions)
    segments(PositionsMb, rep(0,nbG), PositionsMb, Values, col=ValuesColors)
    if(PlottingGenes) {
        segments(PositionsMb, rep(-0.1,nbG), PositionsMb, rep(-1,nbG), col=GeneColors)
    }

}

### Example
ChrLength <- 100000000 # Chromosome of 100 Mb
Positions <- sample(1:100000000, 1000, replace=F)
Values <- sample(1:100 / 10, 1000, replace=T)

# Simple Call
ChrDensity(ChrLength, Positions, Values, main="Test plot")


# Plotting genes and adding specific colors for subset of genes

Colors <- rep("black", length(Positions))
Colors[order(Positions)[200:400]] <- "red"
ChrDensity(ChrLength, Positions, Values, PlottingGenes=T, ValuesColors=Colors, GeneColors=Colors, main="Test plot", ylab="Nb of SNPs")