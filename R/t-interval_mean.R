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

two_sample_pooled_t_interval <- function(sample_1, sample_2, conf_level=0.95){
  
  sample_1_mean <- mean(sample_1)
  sample_2_mean <- mean(sample_2)
  s_1 <- sd(sample_1)
  s_2 <- sd(sample_2)
  n <- length(sample_1)
  m <- length(sample_2)
  alpha <- 1 - conf_level
  
  pooled_sample_variance <- ((n-1)*s_1^2+(m-1)*s_2^2)/(n+m-2)
  
  t_quantile <- qt(alpha/2, n+m-2, lower.tail=FALSE)
  margin_error <- t_quantile*sqrt(pooled_sample_variance*(1/n+1/m))
  mean_diff <- sample_1_mean - sample_2_mean
  interval <- c(lower=mean_diff-margin_error, upper=mean_diff+margin_error, margin_error=margin_error)
}