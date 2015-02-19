######## Auxillary code to generate plots ###############
setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/figures/manhattan")

#Code to plot manhattanplot
chrLens    <- read.table("chrLengths_mm10.txt")$V2
chrLens    <- c(0, cumsum(chrLens*1.0))

labPos     <- read.table("chrLengths_mm10.txt")$V2/2.0
labPos     <- chrLens[1:19]+labPos



traits <- c("cpp.diff", 
            "cpp8.1", 
            #"cpp8.t",
            #"act2.t", 
            "act4.t", 
            #"act5.t",
            #"act8.t", 
            "sens", 
            "startle", 
            #"wild", 
            "ppi3", 
            "ppi6", 
            #"ppi12", 
            "habituation", 
            "glucose")

titles <- c("Change in preference for meth paired side", 
            "Time spent on meth paired side on day 8 (5 min)", 
            #"Time spent on meth paired side on day 8 (30 min)", 
            #"Day 2 activity (1 mg/kg meth)", 
            "Day 4 activity (1 mg/kg meth)", 
            #"Day 5 activity (saline)",
            #"Day 8 activity (saline)",
            "Locomotor sensitization to meth", 
            "Acoustic startle response", 
            #"Wildness (number of escapes during cpp)", 
            "Prepulse inhibition (73 dB)",
            "Prepulse inhibition (76 dB)", 
            #"Prepulse inhibition (82 dB)", 
            "Habituation to acoustic startle",
            "Blood glucose levels after a 4h fast")

outfiles <- paste0(traits, ".manhattan.png")


for (index in 1:9) {
        print(paste("Starting trait", traits[index]))
        plotManhattan(traits[index], outfile=outfiles[index], title=titles[index])
        print(paste("Done with trait", traits[index]))
}







plotManhattan <- function(trait, outfile=NULL, 
                          oddcolor="#ff6040", evencolor="#aab0be", 
                          signifLine=6, signifColor="#000000", title=NULL) {
        if (is.null(outfile)) {
                outfile <- paste0(trait, ".manhattan.png")
        }
        if (is.null(title)) {
                title = trait
        }
        png(file=outfile, height=600, width=1400, bg="transparent")
        pvals <- c()
        positions <- c()
        cols <- c()
        
        for (chrom in 1:19) {
                filename <- paste0(trait, ".chr", chrom, ".assoc.txt")
                print(paste("Starting chromosome", chrom))
                output <- read.table(filename, header=T, as.is=TRUE)
                output$chr <- as.numeric(sub("chr", "", output$chr))
                pvals <- c(pvals, -log10(output$p_lrt))
                positions <- c(positions, chrLens[output$chr]+output$ps)
                if (chrom %% 2) {
                        cols <- c(cols, rep(oddcolor, nrow(output)))
                } else {
                        cols <- c(cols, rep(evencolor, nrow(output)))
                }
                print(paste("Done with chromosome", chrom))
        }
        plot(positions, pvals, col=cols, pch=19, cex=0.78, main=title, ylab='-log10(p-value)', xlab='Position', axes=F, frame.plot=T, cex.lab=1.7, cex.main=1.8, font.main=1)
        axis(2, cex.axis=1.5)
        axis(1, at=labPos, labels=1:19, tick=F, cex.axis=1.5)
        abline(h=signifLine, col=signifColor, lty=2)
        dev.off()
}




# plotZoom <- function(trait, chrom, start, end, outfile=NULL, dotColor="firebrick", signifLine=6, signifColor="darkorange1", title=NULL, scaling=1e6, scalename="(Mb)") {
#     if (is.null(outfile)) {
#         outfile <- paste0(trait, ".chr", chrom, ".", as.character(start), ".", as.character(end), ".zoom.jpg")
#     }
#     if (is.null(title)) {
#         title = trait
#     }
#     jpeg(file=outfile, quality=95)
#     filename <- paste0("/home/shyamg/projects/Mouse/AIL/LgSm-DataProcessing/output/", trait, ".chr", chrom, ".assoc.txt")
#     print(paste("Starting chromosome", chrom))
#     output <- read.table(filename, header=T)
#     pvals <- -log10(output$p_lrt)
#     positions <- output$ps/scaling
#     print(paste("Done with chromosome", chrom))
#     ### Code to plot
#     start = start/scaling
#     end = end/scaling
#     plot(positions, pvals, col=dotColor, pch=19, cex=0.75, main=title, ylab='-log10(p-value)', xlab=paste('Position', scalename), axes=F, frame.plot=T, cex.lab=1.1, xlim=c(start, end))
#     axis(2, cex.axis=1.2)
#     axis(1, cex.axis=1.2)
#     abline(h=signifLine, col=signifColor, lty=2)
#     dev.off()        
# }






