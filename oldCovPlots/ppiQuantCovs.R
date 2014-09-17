panels <-
      list( ppi.sa= list (pheno="avg_ppi",   cov="sire.age",    col="firebrick"),
            ppi.w= list (pheno="avg_ppi",   cov="ppi.weight",  col="dodgerblue"),
            ppi.a= list(pheno="avg_ppi",    cov="ppi.age",     col="darkblue"),
            ppi.cpp= list (pheno="avg_ppi",   cov="cpp.diff",    col="green"),
            ppi.sz= list (pheno="avg_ppi",   cov="sens",        col="lightgreen"),
            ppi.s= list (pheno="avg_ppi",   cov="startle",        col="orange"),
            
            startle.sa= list (pheno="startle",   cov="sire.age",    col="firebrick"),
            startle.a= list(pheno="startle",    cov="ppi.age",     col="darkblue"),
            startle.w= list (pheno="startle",   cov="ppi.weight",  col="dodgerblue"),
            startle.cpp= list (pheno="startle",   cov="cpp.diff",    col="green"),
            startle.sz= list (pheno="startle",   cov="sens",        col="lightgreen"),
            
            hab.sa= list (pheno="habituation",   cov="sire.age",    col="firebrick"),
            hab.a= list(pheno="habituation",    cov="ppi.age",     col="darkblue"),
            hab.w= list (pheno="habituation",   cov="ppi.weight",  col="dodgerblue"),
            hab.cpp= list (pheno="habituation",   cov="cpp.diff",    col="green"),
            hab.sz= list (pheno="habituation",   cov="sens",        col="lightgreen"))


phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno      <- ppidata
n     <- nrow(pheno)
pheno <- transform(pheno,ppi.age = ppi.age + rnorm(n,sd = 0.5))
pheno <- transform(pheno,sire.age = sire.age + rnorm(n,sd = 0.5))



for (panel in names(panels)) {      
      r         <- panels[[panel]]
      phenotype <- r$pheno
      covariate <- r$cov
      panel.col <- r$col
      data        <- pheno[c(phenotype,covariate)]
      #data <- cbind(data, pheno)
      names(data) <- c("y","x")
      model <- lm(y ~ x,data)
      pve   <- summary(model)$r.squared
      
      filename=paste0(covariate,"_", phenotype, ".png", sep="")
      png(file=filename, height=300,width=300)
      
      print(ggplot(data, aes(x=x,y=y)) +
                  geom_point(color=panel.col)+
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

