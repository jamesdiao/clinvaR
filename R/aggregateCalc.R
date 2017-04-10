#' Aggregate Allele Frequencies Using Independence Assumption
#'
#' This function computes an aggregated allele frequency: P(any pathogenic allele) across 
#' a given dataset, at the given locations, with the given inheritance patterns, across
#' ancestral groups. 
#' @usage aggregateCalc(input, superpop, item, dataset, loc, inherit)
#' @examples aggregateCalc(input = merged_1000g, superpop = 'AFR', item = 'BRCA2', 
#' dataset = '1000 Genomes', loc = rep(T, nrow(merged_1000g)), inherit.use = inheritance.gene))

aggregateCalc <- function(input, superpop, item, dataset, loc, inherit) {
  find = sprintf("AF_%s",toupper(dataset))
  if (superpop!=dataset) 
    find = paste(find, superpop, sep = "_")
  # Aggregation by calculation + ind assumption
  freq <- input[loc,find] %>% unlist %>% as.numeric #vector of all allele frequencies
  if (inherit[item] == "AR")
    freq <- freq^2       #AR: prob of having a non-reference site at BOTH alleles
  if (inherit[item] %in% c("AD","SD"))
    freq <- 1-(1-freq)^2 #AD/SD: prob of having a non-reference site at EITHER allele
  if (inherit[item] == "XL")
    freq <- freq         #XL: prob of having a non-reference site at ONE allele (male)
  freq <- 1-prod(1-freq[!is.na(freq)]) # prob of having a non-reference site in 1 chrom
  freq <- freq_CI(freq, 2*sample_size[[paste0("Cohort_",dataset)]][superpop], 0.95)
  return(freq)
}