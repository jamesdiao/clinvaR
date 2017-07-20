#' Compute One-Tailed Beta Distribution Confidence Interval from Observed Alleles
#'
#' This function computes a one-Tailed Beta distribution confidence interval from observed alleles
#' @usage freq_CI(freq, n, cutoff)
#' @param freq numeric; observed frequency
#' @param cutoff numeric; upper tailed confidence range. Default is 95% = 0.95. 
#' @param n numeric; number of observations, equal to the number of observed alleles 
#' (twice the number of observed individuals). This is hard coded using dataset FAQs for 
#' the aggregateCount() and aggregateCalc() functions
#' @examples freq_CI(freq = 0.05, n = 2*2504, cutoff = 0.95)
#' @export

freq_CI <- function(freq, n, cutoff) {
  if (missing(cutoff)) 
    cutoff <- 0.95
  obs_success = freq * n
  obs_failure = n-obs_success
  qbeta(p = cutoff, shape1 = obs_success+1, shape2 = obs_failure+1, lower.tail = F) %>% 
    return()
}
