############ binary.from.categorical ###########################
# This function by Peter Carbonetto turns categorical variables
# (usually covariates) into separate binary phenotypes. It
# returns a df with 1 col for each factor level.
################################################################

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






