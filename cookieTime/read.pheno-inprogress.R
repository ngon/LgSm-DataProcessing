

# Create binary covariates
# Returns a data frame with one column for each level of factor x. 
binary.from.categorical <- function (x, col.names = NULL) {
  # Create a binary factor for each value of the categorical variable.
  d <- list()
  for (value in levels(x))
    d[[value]] <- as.integer(x == value)
  
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

source("../cookieTime/r")
phenotypes <- read.pheno("allData.txt")

# make binary phenotypes from some factor variables
binaryPhenos <- cbind(binary.from.categorical(phenotypes$gen, 
                                       col.names=paste0("is.gen", 50:56)),
               binary.from.categorical(phenotypes$coat, 
                                       col.names=paste0("is.coat", c("A","B","W"))),
               binary.from.categorical(phenotypes$batch, 
                                       col.names=paste0("is.batch", 1:22)),
               binary.from.categorical(phenotypes$cpp.box, 
                                       col.names=paste0("is.cpp.box",1:12)),
               binary.from.categorical(phenotypes$ppi.box, 
                                       col.names=paste0("is.ppi.box",1:5))
               )

# combine with phenotypes data frame
phenotypes <- cbind(phenotypes, binaryPhenos)

# write to an all-in-one file
write.table(phenotypes, "allDataWithBinaryPhenotypes.txt", sep="\t", row.names=F,
            col.names=T)

# make covariate files
one <- data.frame(rep(1, length(phenotypes$id)))
covariates <- cbind(one, phenotypes[1:36], phenotypes[110:116], phenotypes[118:119],
                    phenotypes[121:122], phenotypes[125:128], phenotypes[153:202])
names(covariates)[1] <- "one"

write.table(names(covariates), "cov.names.txt", sep="\t", row.names=F)
write.table(covariates, "cov.noHeader.txt", sep="\t", row.names=F, col.names=F)

# make pheno files
traits <- cbind(phenotypes[c(1, 25:109, 117, 120, 122:124, 126, 129:153, 161:163)])
write.table(names(traits), "pheno.names.txt", sep="\t", row.names=F)
write.table(traits, "pheno.noHeader.txt", sep="\t", row.names=F, col.names=F)




