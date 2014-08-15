

      panels <-
            list(act1.1    = list(pheno="act1.1",    cov="sire.age",    col="firebrick"),
                 act1.t    = list(pheno="act1.t",    cov="sire.age",    col="firebrick"),
                 act2.1    = list(pheno="act2.1",    cov="sire.age",    col="firebrick"),
                 act2.t    = list(pheno="act2.t",    cov="sire.age",    col="firebrick"),
                 act3.1    = list(pheno="act3.1",    cov="sire.age",    col="firebrick"),
                 act3.t    = list(pheno="act3.t",    cov="sire.age",    col="firebrick"),
                 act4.1    = list(pheno="act4.1",    cov="sire.age",    col="firebrick"),
                 act4.t    = list(pheno="act4.t",    cov="sire.age",    col="firebrick"),
                 act5.1    = list(pheno="act5.1",    cov="sire.age",    col="firebrick"),
                 act5.t    = list(pheno="act5.t",    cov="sire.age",    col="firebrick"),
                 act8.1    = list(pheno="act8.1",    cov="sire.age",    col="firebrick"),
                 act8.t    = list(pheno="act8.t",    cov="sire.age",    col="firebrick"),
                 
                 act1.1c    = list(pheno="act1.1",    cov="d1.weight",    col="dodgerblue"),
                 act1.tc    = list(pheno="act1.t",    cov="d1.weight",    col="dodgerblue"),
                 act2.1c    = list(pheno="act2.1",    cov="d2.weight",    col="dodgerblue"),
                 act2.tc    = list(pheno="act2.t",    cov="d2.weight",    col="dodgerblue"),
                 act3.1c    = list(pheno="act3.1",    cov="d3.weight",    col="dodgerblue"),
                 act3.tc    = list(pheno="act3.t",    cov="d3.weight",    col="dodgerblue"),
                 act4.1c    = list(pheno="act4.1",    cov="d4.weight",    col="dodgerblue"),
                 act4.tc    = list(pheno="act4.t",    cov="d4.weight",    col="dodgerblue"),
                 act5.1c    = list(pheno="act5.1",    cov="d5.weight",    col="dodgerblue"),
                 act5.tc    = list(pheno="act5.t",    cov="d5.weight",    col="dodgerblue"),
                 act8.1c    = list(pheno="act8.1",    cov="d8.weight",    col="dodgerblue"),
                 act8.tc    = list(pheno="act8.t",    cov="d8.weight",    col="dodgerblue"),
                 
                 act1.1d    = list(pheno="act1.1",    cov="cpp.age",    col="darkblue"),
                 act1.td    = list(pheno="act1.t",    cov="cpp.age",    col="darkblue"),
                 act2.1d    = list(pheno="act2.1",    cov="cpp.age",    col="darkblue"),
                 act2.td    = list(pheno="act2.t",    cov="cpp.age",    col="darkblue"),
                 act3.1d    = list(pheno="act3.1",    cov="cpp.age",    col="darkblue"),
                 act3.td    = list(pheno="act3.t",    cov="cpp.age",    col="darkblue"),
                 act4.1d    = list(pheno="act4.1",    cov="cpp.age",    col="darkblue"),
                 act4.td    = list(pheno="act4.t",    cov="cpp.age",    col="darkblue"),
                 act5.1d    = list(pheno="act5.1",    cov="cpp.age",    col="darkblue"),
                 act5.td    = list(pheno="act5.t",    cov="cpp.age",    col="darkblue"),
                 act8.1d    = list(pheno="act8.1",    cov="cpp.age",    col="darkblue"),
                 act8.td    = list(pheno="act8.t",    cov="cpp.age",    col="darkblue"))
             


      

phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno      <- natalia
pn     <- nrow(pheno)
pheno <- transform(pheno,cpp.age = cpp.age + rnorm(n,sd = 0.5))



for (panel in names(panels)) {      
      r         <- panels[[panel]]
      phenotype <- r$pheno
      covariate <- r$cov
      #panel.col <- r$col
      data        <- pheno[c(phenotype,covariate)]
      data <- cbind(data, pheno$sex)
      names(data) <- c("y","x", "sex")
      model <- lm(y ~ x + sex,data)
      pve   <- summary(model)$adj.r.squared
      
      filename=paste0("sex_",covariate,"_", phenotype, ".png", sep="")
      png(file=filename, height=300,width=300)
      
      print(ggplot(pheno, aes(x,y, color=box)) +
            geom_point()+
            geom_smooth(method=lm, se=TRUE) +
            xlab(paste0(covariate, " (PVE= ", round(100*pve, digits=2), "%)"))+
            ylab(paste0(phenotype))+
            guides(color=FALSE)+
            theme(axis.title.x = element_text(colour="black", size=16),
                  axis.text = element_text(colour="black", size=12),
                  axis.title.y = element_text(colour="black", size=16),
                  plot.title = element_text(colour="black", size=12),
                  legend.position="none")) 
      
      dev.off()      
}
 
