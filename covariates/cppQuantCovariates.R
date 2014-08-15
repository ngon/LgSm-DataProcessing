

panels <-
      list(cpp1.1    = list(pheno="cpp1.1",    cov="sire.age",    col="firebrick"),
           cpp1.t    = list(pheno="cpp1.t",    cov="sire.age",    col="firebrick"),
           cpp8.1    = list(pheno="cpp8.1",    cov="sire.age",    col="firebrick"),
           cpp8.t    = list(pheno="cpp8.t",    cov="sire.age",    col="firebrick"),
           
 
           cpp1.1c    = list(pheno="cpp1.1",    cov="d1.weight",    col="dodgerblue"),
           cpp1.tc    = list(pheno="cpp1.t",    cov="d1.weight",    col="dodgerblue"),
           cpp8.1c    = list(pheno="cpp8.1",    cov="d8.weight",    col="dodgerblue"),
           cpp8.tc    = list(pheno="cpp8.t",    cov="d8.weight",    col="dodgerblue"),
           
           cpp1.1d    = list(pheno="cpp1.1",    cov="cpp.age",    col="darkblue"),
           cpp1.td    = list(pheno="cpp1.t",    cov="cpp.age",    col="darkblue"),
           cpp8.1d    = list(pheno="cpp8.1",    cov="cpp.age",    col="darkblue"),
           cpp8.td    = list(pheno="cpp8.t",    cov="cpp.age",    col="darkblue"),
           
                    
           
           )





phenotypes <- unique(sapply(panels,function(x)x$pheno))
pheno      <- natalia
pn     <- nrow(pheno)
pheno <- transform(pheno,cpp.age = cpp.age + rnorm(n,sd = 0.5))


effectSizes<- function(phenotype, covariate) {
      
}
for (panel in names(panels)) {      
      r         <- panels[[panel]]
      phenotype <- r$pheno
      covariate <- r$cov
      panel.col <- r$col
      data        <- pheno[c(phenotype,covariate)]
      names(data) <- c("y","x")
      model <- lm(y ~ x,data)
      pve   <- summary(model)$r.squared
      
      filename=paste0(covariate,"_", phenotype, ".png", sep="")
      png(file=filename, height=300,width=300)
      
      print(ggplot(data, aes(x,y)) +
                  geom_point(color=panel.col)+
                  geom_smooth(method=lm, se=TRUE) +
                  xlab(paste0(covariate, " (PVE= ", round(100*pve, digits=2), "%"))+
                  ylab(paste0(phenotype))+
                  theme(axis.title.x = element_text(colour="black", size=16),
                        axis.text = element_text(colour="black", size=12),
                        axis.title.y = element_text(colour="black", size=16),
                        plot.title = element_text(colour="black", size=12))) 
      
      dev.off()      
}

