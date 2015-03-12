

######### FUNCTIONS ################

read.pheno <- function(file){
# read in tab deliminted phenotype file  
pheno <- read.table(file, sep="\t", header=T, stringsAsFactors=F)

# convert some columns to factors
pheno <- transform(pheno,
                   id       = as.character(id),
                   gen      = factor(gen),
                   cc       = factor(cc),
                   wean.age = factor(wean.age),
                   sex      = factor(sex),
                   batch    = factor(batch),
                   cpp.box  = factor(cpp.box),
                   ppi.box  = factor(ppi.box),
                   comerr1  = factor(comerr1),
                   comerr8  = factor(comerr8),
                   brain    = factor(brain)
                   )

# cc has information about the coat color and ear tag. 
# extract information about coat color only; make new variable coat.
pheno$coat <- pheno$cc
levels(pheno$coat)[levels(pheno$coat)=="AGL"] <- "A"
levels(pheno$coat)[levels(pheno$coat)=="AGR"] <- "A"
levels(pheno$coat)[levels(pheno$coat)=="AGRL"] <- "A"
levels(pheno$coat)[levels(pheno$coat)=="AL"] <- "A"
levels(pheno$coat)[levels(pheno$coat)=="AR"] <- "A"
levels(pheno$coat)[levels(pheno$coat)=="ARL"] <- "A"
levels(pheno$coat)[levels(pheno$coat)=="WL"] <- "W"
levels(pheno$coat)[levels(pheno$coat)=="WR"] <- "W"
levels(pheno$coat)[levels(pheno$coat)=="WRL"] <- "W"
levels(pheno$coat)[levels(pheno$coat)=="BL"] <- "B"
levels(pheno$coat)[levels(pheno$coat)=="BR"] <- "B"
levels(pheno$coat)[levels(pheno$coat)=="BRL"] <- "B"


# convert integer columns to numeric columns
pheno <- transform(pheno,
                   act1.1   = as.numeric(act1.1),
                   act1.2   = as.numeric(act1.2),
                   act1.3   = as.numeric(act1.3),
                   act1.4   = as.numeric(act1.4),
                   act1.5   = as.numeric(act1.5),
                   act1.6   = as.numeric(act1.6),
                   act1.t   = as.numeric(act1.t),
                   act2.1   = as.numeric(act2.1),
                   act2.2   = as.numeric(act2.2),
                   act2.3   = as.numeric(act2.3),
                   act2.4   = as.numeric(act2.4),
                   act2.5   = as.numeric(act2.5),
                   act2.6   = as.numeric(act2.6),
                   act2.t   = as.numeric(act2.t),
                   act3.1   = as.numeric(act3.1),
                   act3.2   = as.numeric(act3.2),
                   act3.3   = as.numeric(act3.3),
                   act3.4   = as.numeric(act3.4),
                   act3.5   = as.numeric(act3.5),
                   act3.6   = as.numeric(act3.6),
                   act3.t   = as.numeric(act3.t),
                   act4.1   = as.numeric(act4.1),
                   act4.2   = as.numeric(act4.2),
                   act4.3   = as.numeric(act4.3),
                   act4.4   = as.numeric(act4.4),
                   act4.5   = as.numeric(act4.5),
                   act4.6   = as.numeric(act4.6),
                   act4.t   = as.numeric(act4.t),
                   act5.1   = as.numeric(act5.1),
                   act5.2   = as.numeric(act5.2),
                   act5.3   = as.numeric(act5.3),
                   act5.4   = as.numeric(act5.4),
                   act5.5   = as.numeric(act5.5),
                   act5.6   = as.numeric(act5.6),
                   act5.t   = as.numeric(act5.t),
                   act8.1   = as.numeric(act8.1),
                   act8.2   = as.numeric(act8.2),
                   act8.3   = as.numeric(act8.3),
                   act8.4   = as.numeric(act8.4),
                   act8.5   = as.numeric(act8.5),
                   act8.6   = as.numeric(act8.6),
                   act8.t   = as.numeric(act8.t),
                   
                   sc1.1    = as.numeric(sc1.1),
                   sc1.2    = as.numeric(sc1.2),
                   sc1.3    = as.numeric(sc1.3),
                   sc1.4    = as.numeric(sc1.4),
                   sc1.5    = as.numeric(sc1.5),
                   sc1.6    = as.numeric(sc1.6),
                   sc1.t    = as.numeric(sc1.t),
                   
                   sc8.1    = as.numeric(sc8.1),
                   sc8.2    = as.numeric(sc8.2),
                   sc8.3    = as.numeric(sc8.3),
                   sc8.4    = as.numeric(sc8.4),
                   sc8.5    = as.numeric(sc8.5),
                   sc8.6    = as.numeric(sc8.6),
                   sc8.t    = as.numeric(sc8.t),
                 
                   sens     = as.numeric(sens),
                   sens1    = as.numeric(sens1),
                   sens2    = as.numeric(sens2),
                   sens3    = as.numeric(sens3),
                   sens4    = as.numeric(sens4),
                   sens5    = as.numeric(sens5),
                   sens6    = as.numeric(sens6),
                   
                   sum.nostim = as.numeric(sum.nostim)
                   
                   )
return(pheno)
}


# Create binary covariates
# Returns a data frame with one column for each level of factor x. 
binary.from.categorical <- function (x, col.names = NULL) {
  # Create a binary factor for each value of the categorical variable.
  d <- list()
  for (value in levels(x))
    d[[value]] <- factor(as.integer(x == value))
  
  # Convert the list to a data frame, and adjust the column names, if
  # requested.
  d <- data.frame(d,check.names = FALSE)
  if (!is.null(col.names))
    names(d) <- col.names
  
  # Output the newly created data frame.
  return(d)
}

########################
setwd("C:/Users/Administrator/Desktop/Scripts/LgSm-DataProcessing/dataFiles")


pheno<- read.pheno("allData.txt")

# make binary phenotypes from some factor variables
pheno <- cbind(pheno,
               binary.from.categorical(pheno$gen, col.names=paste0("gen",50:56)),
               binary.from.categorical(pheno$batch, col.names=paste0("batch",1:22)),
               binary.from.categorical(pheno$cpp.box, col.names=paste0("cpp.box",1:12)),
               binary.from.categorical(pheno$ppi.box, col.names=paste0("ppi.box",1:5)),
               binary.from.categorical(pheno$batch, col.names=paste0("coat",1:3))
               )

# not sure why, but running the command above creates 19 NA variables
# at the end of the data frame pheno. i'm deleting them. 
pheno <- pheno[,-c(202:220)]



















# This function prepares the phenotype data for QTL mapping.
prepare.pheno <- function (pheno) {
   
  # The deaf mice add "noise" to the PPI phenotypes. remove these samples 
  # from the analysis. My cutoff is 10, which is what Peter used for CFWs.
  ppi.phenotypes <-
    c("ppi3", "ppi6", "ppi12", "ppi3.logit", "ppi6.logit", "ppi12.logit", 
      "startle", "habituation")
  rows <- which(pheno$startle < 10) 
  pheno[rows,ppi.phenotypes] <- NA
  
   # Return the updated phenotype data table.
  return(pheno)
}


