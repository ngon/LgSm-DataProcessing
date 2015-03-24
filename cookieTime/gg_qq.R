
############################### gg_qq ############################## 
#
#This function takes the residuals of a linear model as its input 
# and draws a qqplot using the package 'ggplot2' which can be used to 
# identify outliers. It returns a plot with confidence intervals and a 
# data frame of residuals that fall outside the confidence intervals. 
#
# Part of this code comes from car:::qqPlot. Source: StackOverflow


# example call
# x <- rstudent(model); gg_qq(x, trait="cpp1.t")

gg_qq <- function(x, trait, distribution = "norm", ..., line.estimate = NULL, 
                  conf = 0.99, labels = names(x)){
           
  # prepare residuals for qqplotting
  q.function <- eval(parse(text = paste0("q", distribution)))
  d.function <- eval(parse(text = paste0("d", distribution)))
  x <- na.omit(x)
  ord <- order(x)
  n <- length(x)
  P <- ppoints(length(x))
  df <- data.frame(ord.x = x[ord], z = q.function(P, ...))
    
  # get qqline
  if(is.null(line.estimate)){
    Q.x <- quantile(df$ord.x, c(0.25, 0.75))
    Q.z <- q.function(c(0.25, 0.75), ...)
    b <- diff(Q.x)/diff(Q.z)
    coef <- c(Q.x[1] - b * Q.z[1], b)
  } else {
    coef <- coef(line.estimate(ord.x ~ z))
  }
  
  ## calculate the confidence intervals
  zz <- qnorm(1 - (1 - conf)/2)
  SE <- (coef[2]/d.function(df$z)) * sqrt(P * (1 - P)/n)
  fit.value <- coef[1] + coef[2] * df$z
  df$upper <- fit.value + zz * SE
  df$lower <- fit.value - zz * SE
  
  ##
  norm.r <- (df$ord.x - mean(df$ord.x))/sqrt(var(df$ord.x))
  oeQ <- data.frame(upper=  mean(-3 < norm.r & norm.r < 3),
                    lower= 1-(mean(-3 < norm.r & norm.r < 3)))
  
    ## label the points outside the confidence intervals with their indexes
  if(!is.null(labels)){ 
    df$ci.out <- ifelse(df$ord.x > df$upper | df$ord.x < df$lower, labels[ord],"")
  }
   df$sd3.out <- ifelse(df$ord.x > quantile(df$ord.x, oeQ$upper) | df$ord.x < quantile(df$ord.x, oeQ$lower), labels[ord], "")

  
  ## qqplot with confidence intervals - ggplot2
  p <- ggplot(df, aes(x=z, y=ord.x)) +
    geom_point() + 
    geom_abline(intercept = coef[1], slope = coef[2]) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha=0.2) +
    labs(title = paste0(trait, "_outliers"),
         x = "Theoretical quantiles",
         y = "Standardized residuals") +
    theme(legend.position="none")
  if(!is.null(labels)) p <- p + geom_text( aes(label = ci.out, size=8, angle=90,
                                               vjust=2, hjust=.85))
  p <- p + geom_text( aes(label = sd3.out, size=8, angle=90,
                          vjust=2, hjust=.85, colour="red"))
  png(filename=paste0("QQ-", trait, ".png"), 
       height =350, width=1200, units="px",pointsize=8)
  print(p)
  coef
  dev.off()  
  
  # one: return indexes for values outside the specified confidence interval
  # two: return indexes for values > 3 sd from the mean
  one <- df[df$ci.out >= 0,]
  two <- df[df$sd3.out >=0,]
  return(as.list(merge(one, two, all=TRUE )))
   
}
  
  