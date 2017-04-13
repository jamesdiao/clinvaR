#' Aggregate Allele Frequencies by Counting in 1000 Genomes
#'
#' This function computes an aggregated allele frequency: P(any pathogenic allele) across 
#' a given dataset, at the given locations, with the given inheritance patterns, across
#' ancestral groups. 
#' 
#' @usage aggregateCount(input, superpop, item, dataset, loc, inherit, sample_size)
#' @examples aggregateCount(input = merged_1000g, superpop = 'AFR', item = 'BRCA2', dataset = '1000 Genomes', 
#' loc = rep(T, nrow(merged_1000g)), inherit.use = inheritance.gene, sample_size = sample_size)
#' @export

aggregateCount <- function(input, superpop, item, dataset, loc, inherit, sample_size) {
  # Aggregation by counting
  front_cols <- 1:(grep("HG00096",colnames(input))-1)
  find <- (1:ncol(input))[-front_cols]
  if (superpop != dataset) 
    find <- length(front_cols)+which(map$super_pop==superpop)
  if (inherit[item] %in% c("AD","SD"))
    reduced_input <- input[loc, find]
  if (inherit[item] == "AR")
    reduced_input <- input[loc, find]-1 #Looking for 2s
  if (inherit[item] == "XL") {
    male <- length(front_cols)+which(map$gender=="male")
    reduced_input <- input[loc,intersect(find,male)]
  }
  freq <- mean(apply(reduced_input, 2, function(col) any(col>=1))) %>% 
    freq_CI(2*sample_size[["Cohort_1000G"]][superpop], 0.95)
  return(freq)
}