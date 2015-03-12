
getP.plotQ <- function(directory=getwd(), trait){
      
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
      jpeg(filename=paste(trait, "qqplot.jpg", sep=""))
      
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




