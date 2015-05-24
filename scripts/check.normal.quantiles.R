check.normal.quantiles <- function (x) {
  
  # Discard the missing values.
  x <- x[!is.na(x)]
  
  # Transform the observations so that they have zero mean and unit
  # variance.
  x <- (x - mean(x))/sqrt(var(x))
  
  # Create a data frame giving the observed and expected proportion of
  # samples within 1, 2 and 3 standard deviations of the mean.
  a<- data.frame(exp = c(pnorm(1) - pnorm(-1),
                            pnorm(2) - pnorm(-2),
                            pnorm(3) - pnorm(-3)),
                    obs = c(mean(-1 < x & x < 1),
                            mean(-2 < x & x < 2),
                            mean(-3 < x & x < 3)),
                    row.names = c("sd1","sd2","sd3"))
   return(a)
}