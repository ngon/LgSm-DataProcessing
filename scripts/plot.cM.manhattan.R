### PURPOSE: plot association results with cM position on x-axis.

plot.cM <- function(trait, outfile=NULL,
                          oddcolor="#ff6040", evencolor="#aab0be",
                          signifLine=5, signifColor="#000000", title=NULL) {
    if (is.null(outfile)) {
        outfile <- paste0(trait, ".cMhattan.pdf")
    }
    if (is.null(title)) {
        title = trait
    }

    for (chrom in 1:19) {
        filename <- paste0("/group/palmer-lab/AIL/qtlmapping/output/cM.",
                           trait, ".chr", chrom, ".assoc.txt")

        print(paste("Starting chromosome ", chrom))
        output <- read.table(filename, header=T, as.is=TRUE, sep="\t")
       # output$chr <- as.numeric(sub("chr", "", output$chr))
        pvals <- (-log10(output$p_lrt))
        cms <- output$cM
        if (chrom %% 2) {
            cols <- oddcolor
        } else {
            cols <- evencolor
        }
        print(paste0("Done with chromosome ", chrom, ": ready to plot."))


       pdf(file=outfile, height=2.25, width=5, bg="transparent")

       plot(cms, pvals, col=cols, pch=19, cex=0.78, main=title, ylab='-log10(p-value)',
            xlab='Genetic distance (cM)', axes=T, frame.plot=T, cex.lab=1.7,
            cex.main=1.8, font.main=1)
       # axis(2, cex.axis=1.5)
      # axis(1, at=labPos, labels=1:19, tick=F, cex.axis=1.5)
       abline(h=signifLine, col=signifColor, lty=2)
       dev.off()

    }

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

#traits <- vector("list", length=96)

traits <- c("ppi3.logit", "ppi6.logit", "ppi12.logit", "habituation", "startle",
                      "sc1.t", "sc1.1", "sc1.2", "sc1.3", "sc1.4", "sc1.5", "sc1.6",
                      "sc8.t", "sc8.1", "sc8.2", "sc8.3", "sc8.4", "sc8.5", "sc8.6",
                      "cpp1.t", "cpp1.1", "cpp1.2", "cpp1.3", "cpp1.4", "cpp1.5", "cpp1.6",
                      "cpp8.t", "cpp8.1", "cpp8.2", "cpp8.3", "cpp8.4", "cpp8.5", "cpp8.6",
                      "cpp.diff","cpp.diff.p","cpp.diff1", "cpp.diff2", "cpp.diff3",
                      "cpp.diff4","cpp.diff5","cpp.diff6",
                      "act1.t", "act1.1", "act1.2", "act1.3", "act1.4", "act1.5", "act1.6",
                      "act2.t", "act2.1", "act2.2", "act2.3", "act2.4", "act2.5", "act2.6",
                      #"act3.t", "act3.1", "act3.2", "act3.3", "act3.4", "act3.5", "act3.6",
                      "act4.t", "act4.1", "act4.2", "act4.3", "act4.4", "act4.5", "act4.6",
                      "act5.t", "act5.1", "act5.2", "act5.3", "act5.4", "act5.5", "act5.6",
                      "act8.t", "act8.1", "act8.2", "act8.3", "act8.4", "act8.5", "act8.6",
                      "sens",  "sens1", "sens2", "sens3", "sens4", "sens5", "sens6",
                      "wild.binary", "tail","glucose",
                      "is.coatA", "is.coatB", "is.coatW")




titles <- c("Prepulse inhibition (73 dB)","Prepulse inhibition (76 dB)",
            "Prepulse inhibition (82 dB)", "Habituation to acoustic startle",
            "Acoustic startle response",

            "Side changes on Day 1", "Side changes on Day 1 (0-5 min)",
            "Side changes on Day 1 (5-10 min)", "Side changes on Day 1 (10-15 min)",
            "Side changes on Day 1 (15-20 min)", "Side changes on Day 1 (20-25 min)",
            "Side changes on Day 1 (25-30 min)",

            "Side changes on Day 8", "Side changes on Day 8 (0-5 min)",
            "Side changes on Day 8 (5-10 min)", "Side changes on Day 8 (10-15 min)",
            "Side changes on Day 8 (15-20 min)", "Side changes on Day 8 (20-25 min)",
            "Side changes on Day 8 (25-30 min)",

            "Initial preference", "Initial preference (0-5 min)",
            "Initial preference (5-10 min)",  "Initial preference (10-15 min)",
            "Initial preference (15-20 min)",  "Initial preference (20-25 min)",
            "Initial preference (25-30 min)",
            "Final preference", "Final preference (0-5 min)",
            "Final preference (5-10 min)",  "Final preference (10-15 min)",
            "Final preference (15-20 min)",  "Final preference (20-25 min)",
            "Final preference (25-30 min)",
            "Change in preference for meth paired side (seconds)",
            "Change in preference for meth paired side (% time)",
            "Change in preference for meth paired side (0-5 min)",
            "Change in preference for meth paired side (5-10 min)",
            "Change in preference for meth paired side (10-15 min)",
            "Change in preference for meth paired side (15-20 min)",
            "Change in preference for meth paired side (20-25 min)",
            "Change in preference for meth paired side (25-30 min)",

            "Day 1 activity (saline)", "Day 1 activity (saline, 0-5 min)",
            "Day 1 activity (saline, 5-10 min)", "Day 1 activity (saline, 10-15 min)",
            "Day 1 activity (saline, 15-20 min)", "Day 1 activity (saline, 20-25 min)",
            "Day 1 activity (saline, 25-30 min)",

            "Day 2 activity (meth)", "Day 2 activity (meth, 0-5 min)",
            "Day 2 activity (meth, 5-10 min)", "Day 2 activity (meth, 10-15 min)",
            "Day 2 activity (meth, 15-20 min)", "Day 2 activity (meth, 20-25 min)",
            "Day 2 activity (meth, 25-30 min)",

           # "Day 3 activity (saline)", "Day 3 activity (saline, 0-5 min)",
            #"Day 3 activity (saline, 5-10 min)", "Day 3 activity (saline, 10-15 min)",
            #"Day 3 activity (saline, 15-20 min)", "Day 3 activity (saline, 20-25 min)",
            #"Day 3 activity (saline, 25-30 min)",

            "Day 4 activity (meth)", "Day 4 activity (meth, 0-5 min)",
            "Day 4 activity (meth, 5-10 min)", "Day 4 activity (meth, 10-15 min)",
            "Day 4 activity (meth, 15-20 min)", "Day 4 activity (meth, 20-25 min)",
            "Day 4 activity (meth, 25-30 min)",

            "Day 5 activity (saline)", "Day 5 activity (saline, 0-5 min)",
            "Day 5 activity (saline, 5-10 min)", "Day 5 activity (saline, 10-15 min)",
            "Day 5 activity (saline, 15-20 min)", "Day 5 activity (saline, 20-25 min)",
            "Day 5 activity (saline, 25-30 min)",

            "Day 8 activity (saline)", "Day 8 activity (saline, 0-5 min)",
            "Day 8 activity (saline, 5-10 min)", "Day 8 activity (saline, 10-15 min)",
            "Day 8 activity (saline, 15-20 min)", "Day 8 activity (saline, 20-25 min)",
            "Day 8 activity (saline, 25-30 min)",

            "Locomotor sensitization to meth", "Locomotor sensitization to meth (0-5 min)",
            "Locomotor sensitization to meth (5-10 min)",
            "Locomotor sensitization to meth (10-15 min)",
            "Locomotor sensitization to meth(15-20 min)",
            "Locomotor sensitization to meth (20-25 min)",
            "Locomotor sensitization to meth (25-30 min",

            "Wildness (number of escapes during cpp)",
            "Tail length (cm)",
            "Blood glucose levels after a 4-hour fast",
            "Agouti coat", "Black coat", "White coat"

            )

outfiles <- paste0("/group/palmer-lab/AIL/LgSm-DataProcessing/figures/manhattan/onlyEmpManhattan/", traits, "SM.pdf")






for (index in 1:89) {
    print(paste("Starting trait", traits[index]))
    plotManhattan(traits[index], outfile=outfiles[index], title=titles[index])
    print(paste("Done with trait", traits[index]))
}
