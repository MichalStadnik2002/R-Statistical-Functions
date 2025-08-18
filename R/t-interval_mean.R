t_interval_mean <- function(sample, conf_level=0.95){
  
  sample_mean <- mean(sample)
  s <- sd(sample)
  n <- length(sample)
  alpha <- 1 - conf_level
  
  t_quantile <- qt(alpha/2, n-1, lower.tail=FALSE)
  margin_error <- t_quantile*s/sqrt(n)
  interval <- c(lower=sample_mean-margin_error, upper=sample_mean+margin_error, margin_error=margin_error)
  
  return(interval)
}