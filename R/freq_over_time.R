#' Population frequencies over time
#'
#' This function computes the probability of having at least 1 variant  
#' in the gene list of interest (assuming independence) for a range of ClinVar files
#' 
#' @usage freq_over_time(vcf)
#' freq_over_time(genes)
#' @param vcf data.frame; clinvaR-processed vcf from import_file_1000g()
#' @param genes character; genes of interest.
#' @examples 
#' freq_over_time(get_genes('lmm_hcm'))
#' freq_over_time(import_file_1000g())
#' @export

freq_over_time <- function(genes, vcf, dates) {
  if (missing(vcf) & missing(genes)) {
    stop('Input either VCF from import_file_1000g() or a gene list.')
  }
  if (missing(vcf)) 
    vcf <- import_file_1000g(genes)
  if (missing(dates))
    all_dates <- get_date_list()
  freqs <- lapply(all_dates, function(date) {
    merged_vcf <- merge_clinvar_1000g(clinvar = get_clinvar(date), vcf = vcf) 
    return(1-prod(1-merged_vcf$AF_1000G))
  }) %>% unlist %>% setNames(all_dates)
  plot(all_dates, freqs, type = 'o', xlab = "Dates", ylab = "Frequency", main = "Aggregate Population Frequency")
  return(freqs)
}


