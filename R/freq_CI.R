#' Compute One-Tailed Beta Distribution Confidence Interval from Observed Alleles
#'
#' This function computes a one-Tailed Beta distribution confidence interval from observed alleles
#' @usage freq_CI(freq, n, cutoff)
#' @examples freq_CI(freq = 0.05, n = 2*sample_size, cutoff = 0.95)
#' @export

freq_CI <- function(freq, n, cutoff) {
  obs_success = freq * n
  obs_failure = n-obs_success
  return(qbeta(p = cutoff, shape1 = obs_success+1, shape2 = obs_failure+1, lower.tail = F))
}