cohens_d_calc <- function(x1, x2) {
  m1 <- mean(x1, na.rm = TRUE)
  m2 <- mean(x2, na.rm = TRUE)
  s1 <- var(x1, na.rm = TRUE)
  s2 <- var(x2, na.rm = TRUE)
  pooled_sd <- sqrt((s1 + s2) / 2)
  d <- (m1 - m2) / pooled_sd
  return(d)
}