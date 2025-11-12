jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  return(intersection / union)
}