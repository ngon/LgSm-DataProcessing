


# Functions to remove outliers

# read in pheno data and make all columns numeric.
data <- read.table("phenotypes.txt", sep="\t", header=T)
numericData <- as.data.frame(sapply(data, as.numeric))

# find outliers +3 SD from the mean
find.high.outlier <- function(data, cutoff=3){
  # calculate standard deviation of each col
  sds <- sapply(data, sd, na.rm=TRUE)
  means <- sapply(data, mean, na.rm=TRUE)
  # find cells with value greater than cutoff * sd
  found <- mapply(function(d,m,s){
    which(d > m +cutoff * s)}, data, means, sds)
  found    
  }

# find outliers -3 SD from the mean
find.low.outlier <- function(data, cutoff=3){
  # calculate standard deviation of each col
  sds <- sapply(data, sd, na.rm=TRUE)
  means <- sapply(data, mean, na.rm=TRUE)
  # find cells with value greater than cutoff * sd
  found <- mapply(function(d,m,s){
    which(d < m -cutoff * s)}, data, means, sds)
  found    
}

# remove outlier values and replace with NA
# remove.outlier <- function(data, outliers) {
#   rmout <- mapply(function(d,o){
#     result <- d
#     result[o] <- NA
#     return(result)
#   }, data, outliers)
#   return(as.data.frame(result))
# }  



  
  
  
  
  
  
  
  
  
  
  
  
  
