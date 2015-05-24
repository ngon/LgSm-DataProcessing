setwd("/group/palmer-lab/AIL/qtlmapping/output")

traits <- c("ppi3.logit", "ppi6.logit", "ppi12.logit", "habituation", "startle",
            "sc1.t", "sc1.1", "sc1.2", "sc1.3", "sc1.4", "sc1.5", "sc1.6",
            "sc8.t", "sc8.1", "sc8.2", "sc8.3", "sc8.4", "sc8.5", "sc8.6",
            "cpp1.t", "cpp1.1", "cpp1.2", "cpp1.3", "cpp1.4", "cpp1.5", "cpp1.6",
            "cpp8.t", "cpp8.1", "cpp8.2", "cpp8.3", "cpp8.4", "cpp8.5", "cpp8.6",
            "cpp.diff","cpp.diff.p","cpp.diff1", "cpp.diff2", "cpp.diff3",
            "cpp.diff4","cpp.diff5","cpp.diff6",
            "act1.t", "act1.1", "act1.2", "act1.3", "act1.4", "act1.5", "act1.6",
            "act2.t", "act2.1", "act2.2", "act2.3", "act2.4", "act2.5", "act2.6",
            #"act3.t",
            "act3.1",
            #"act3.2", "act3.3", "act3.4", "act3.5", "act3.6",
            "act4.t", "act4.1", "act4.2", "act4.3", "act4.4", "act4.5", "act4.6",
            "act5.t", "act5.1", "act5.2", "act5.3", "act5.4", "act5.5", "act5.6",
            "act8.t", "act8.1", "act8.2", "act8.3", "act8.4", "act8.5", "act8.6",
            "sens",  "sens1", "sens2", "sens3", "sens4", "sens5", "sens6",
            "wild.binary", "tail","glucose",
            "is.coatA", "is.coatB", "is.coatW")



getP.plotQ <- function(trait){

    directory=getwd()
    # match file names to trait
    fileNames <- vector(length=19)
    chromosomes <- 1:19

    for (i in seq_along(chromosomes)){
        fileNames[i] <- paste(trait,".chr", chromosomes[i], ".assoc.txt", sep="")
    }

    # read only p_lrt columns from selected files and compile them into a single data frame
    print("I'll fetch the files at once, master.")
    chosenFiles <- lapply(file.path(directory, fileNames), read.table, sep="\t",
                          header=T, colClasses = c(rep("NULL", 8), "numeric"))
    pValues <- do.call(what=rbind.data.frame, args=chosenFiles)

    # take the -log10 and get expected and observed values
    pValues <- pValues$p_lrt
    obs = -log10(sort(pValues,decreasing=F))
    exp = -log10( 1:length(obs)/length(obs) )

    print("I'm plotting as fast as I can!")
    jpeg(filename=paste0("/group/palmer-lab/AIL/LgSm-DataProcessing/figures/", trait, ".qq.jpg"))

    p<- plot(exp,obs,pch=19,cex=0.25,
             main=paste(trait, "QQ PLOT", sep=" "),
             xlab=expression(Expected~~-log[10](italic(p))),
             ylab=expression(Observed~~-log[10](italic(p))),
             xlim=c(0,max(exp)),
             ylim=c(0,max(exp)));
    lines(exp,exp,col="red")
    p

    graphics.off()

}

lapply(traits, getP.plotQ)



