median_cl_boot <- function(x, conf = 0.95) {
  x <- x[!is.na(x)]
  n <- length(x)
  # If sample is too small, just return quantiles
  if (n < 10) {
    q <- quantile(x, probs = c(0.5, (1-conf)/2, 1-(1-conf)/2))
    return(data.frame(y = q[1], ymin = q[2], ymax = q[3]))
  }
  # Bootstrap the median
  bmedian <- function(x, i) median(x[i])
  boot.out <- boot(data = x, statistic = bmedian, R = 1000)
  ci <- boot.ci(boot.out, type = "perc", conf = conf)$percent[4:5]
  return(data.frame(y = median(x), ymin = ci[1], ymax = ci[2]))
}