#' Aggregate Allele Frequencies by Counting in 1000 Genomes
#'
#' This function computes an aggregated allele frequency: P(any pathogenic allele) across 
#' a given dataset, at the given locations, with the given inheritance patterns, across ancestral groups. 
#' 
#' @usage aggregate_count(input, superpop, item, dataset, loc, inherit)
#' @param input data frame; merged ClinVar-sequencing dataset. Defaults to merge_clinvar_[dataset](). 
#' @param superpop character; choose from 'AFR','AMR','EAS','EUR','Total'. Defaults to 'Total'. 
#' @param item character; can be an item from GENE, MIM, or MEDGEN. 
#' @param dataset character; choose from '1000 Genomes','ExAC', or 'gnomAD'. 
#' Not case-sensitive. Defaults to 'gnomAD'. 
#' @param loc numeric; subsets variants related to a certain disease-condition. Defaults to rep(T, nrow(input)).  
#' @param inherit vector of named characters; from 'AD', 'SD', AR', 'XL'. 
#' Names are a list of items from GENE, MIM, or MEDGEN. 
#' @examples aggregate_count(input = merged_1000g, superpop = 'AFR', item = 'BRCA2', dataset = '1000 Genomes', 
#' loc = rep(T, nrow(merged_1000g)), inherit.use = inheritance.gene)
#' @export

aggregate_count <- function(input, superpop, item, dataset, loc, inherit) {
  if (missing(dataset)) dataset <- 'gnomAD'
  if (missing(input)) input <- eval(parse(text=sprintf('merge_clinvar_%s()', tolower(dataset)) ))
  if (missing(superpop)) superpop <- 'Total'
  if (missing(loc)) loc <- rep(T, nrow(input))
  
  sample_size <- list()
  sample_size$Cohort_1000G <- table(map$super_pop) %>% 
    as.numeric %>% setNames(super.levels) %>% c("1000G"=2504)
  sample_size$Cohort_EXAC <- c("AFR"=5203,"AMR"=5789,"EAS"=4327,
                               "EUR"=3307+33370,"SAS"=8256,"EXAC"=60252)
  sample_size$Cohort_GNOMAD <- c("AFR"=12942,"AMR"=18237,"EAS"=9472,
                                 "EUR"=5081+13046+63416,"SAS"=15450,"GNOMAD"=137644)
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