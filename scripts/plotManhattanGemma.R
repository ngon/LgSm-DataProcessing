######## Auxillary code to generate plots ###############
#Code to plot manhattanplot
chrLens    <- read.table('/home/shyamg/projects/Mouse/CFW/maps/chrLengths_mm10.txt')$V2
chrLens    <- c(0, cumsum(chrLens*1.0))
labPos     <- (read.table('/home/shyamg/projects/Mouse/CFW/maps/chrLengths_mm10.txt')$V2)/2.0
labPos     <- chrLens[1:19]+labPos

plotManhattan <- function(trait, outfile=NULL, oddcolor="deepskyblue", evencolor="darkgreen", signifLine=6, signifColor="darkorange1", title=NULL) {
    if (is.null(outfile)) {
        outfile <- paste0(trait, ".manhattan.jpg")
    }
    if (is.null(title)) {
        title = trait
    }
    jpeg(file=outfile, quality=95, height=800, width=1600)
    pvals <- c()
    positions <- c()
    cols <- c()
    for (chrom in 1:19) {
        filename <- paste0("/home/shyamg/projects/Mouse/AIL/LgSm-DataProcessing/output/", trait, ".chr", chrom, ".assoc.txt")
        print(paste("Starting chromosome", chrom))
        output <- read.table(filename, header=T)
        pvals <- c(pvals, -log10(output$p_lrt))
        positions <- c(positions, chrLens[output$chr]+output$ps)
        if (chrom %% 2) {
            cols <- c(cols, rep(oddcolor, nrow(output)))
        } else {
            cols <- c(cols, rep(evencolor, nrow(output)))
        }
        print(paste("Done with chromosome", chrom))
    }
    plot(positions, pvals, col=cols, pch=19, cex=0.75, main=title, ylab='-log10(p-value)', xlab='Position', axes=F, frame.plot=T, cex.lab=1.1)
    axis(2, cex.axis=1.2)
    axis(1, at=labPos, labels=1:19, tick=F, cex.axis=1.2)
    abline(h=signifLine, col=signifColor, lty=2)
    dev.off()
}

plotZoom <- function(trait, chrom, start, end, outfile=NULL, dotColor="firebrick", signifLine=6, signifColor="darkorange1", title=NULL, scaling=1e6, scalename="(Mb)") {
    if (is.null(outfile)) {
        outfile <- paste0(trait, ".chr", chrom, ".", as.character(start), ".", as.character(end), ".zoom.jpg")
    }
    if (is.null(title)) {
        title = trait
    }
    jpeg(file=outfile, quality=95)
    filename <- paste0("/home/shyamg/projects/Mouse/AIL/LgSm-DataProcessing/output/", trait, ".chr", chrom, ".assoc.txt")
    print(paste("Starting chromosome", chrom))
    output <- read.table(filename, header=T)
    pvals <- -log10(output$p_lrt)
    positions <- output$ps/scaling
    print(paste("Done with chromosome", chrom))
    ### Code to plot
    start = start/scaling
    end = end/scaling
    plot(positions, pvals, col=dotColor, pch=19, cex=0.75, main=title, ylab='-log10(p-value)', xlab=paste('Position', scalename), axes=F, frame.plot=T, cex.lab=1.1, xlim=c(start, end))
    axis(2, cex.axis=1.2)
    axis(1, cex.axis=1.2)
    abline(h=signifLine, col=signifColor, lty=2)
    dev.off()        
}


traits <- c("cpp8.t", "sens", "act1.t", "act2.t", "act3.t", "act4.t", "act5.t", "act8.t", "avg.ppi", "startle", "glucose", "wild")
titles <- c("Total time spent on meth paired side\non day 8", "Total activity on day 1", "Total activity on day 2", "Total activity on day 3", "Total activity on day 4", "Total activity on day 5", "Total activity on day 8", "Average prepulse inhibition", "Startle response", "Blood glucose", "Wildness")
outfiles <- paste0("/home/shyamg/projects/Mouse/AIL/LgSm-DataProcessing/figures/", traits, ".manhattan.jpg")

for (index in 1:12) {
    print(paste("Starting trait", traits[index]))
    plotManhattan(traits[index], outfile=outfiles[index], title=titles[index])
    print(paste("Done with trait", traits[index]))
}
