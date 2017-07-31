#' Annotate 1000 Genomes Variant VCF using ClinVar
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies. 
#' 
#' @usage annotate_1000g(vcf)
#' annotate_1000g(genes)
#' annotate_1000g(vcf, conflicts = FALSE)
#' @param clinvar data.frame; clinvaR-processed vcf containing ClinVar data. 
#' Defaults to get_clinvar().  
#' @param vcf data.frame; clinvaR-processed vcf containing 1000 genomes sequencing data. 
#' Defaults to importing data from all downloads: vcf <- import_file_1000g(). 
#' @examples 
#' annotate_1000g(genes = get_genes())
#' annotate_1000g(clinvar = get_clinvar(), vcf = import_file_1000g())
#' @export

annotate_1000g <- function(vcf, genes, clinvar, conflicts) {
  if (missing(clinvar)) {
    clinvar <- get_clinvar()
  }
  if (missing(vcf)) {
    if (missing(genes)) {
      vcf <- import_file_1000g()
    } else {
      download_1000g(genes)
      vcf <- import_file_1000g(genes)
    }
  }
  if (missing(conflicts))
    conflicts <- TRUE
  if (conflicts) {
    keep <- clinvar$pathogenic_incl_conflicts
  } else {
    keep <- clinvar$pathogenic_no_conflicts
  }
  inter <- intersect(clinvar$VAR_ID[keep], vcf$VAR_ID)
  clinvar_merged <- clinvar[(clinvar$VAR_ID %in% inter),] %>% arrange(VAR_ID) %>% select(-GENE)
  vcf_merged <- vcf[vcf$VAR_ID %in% inter,] %>% arrange(VAR_ID)
  front_cols <- 1:(grep("HG00096",colnames(vcf))-1)
  vcf_merged <- data.frame(vcf_merged[,c("GENE","AF_1000G")], 
        clinvar_merged,vcf_merged[,-front_cols])
  class(vcf_merged) <- c("annotated_vcf", "plain_vcf", class(vcf_merged))
  return(vcf_merged)
}

#' Compute and Plot Aggregate Allele Frequencies for 1000 Genomes
#'
#' This function takes an annotated ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of counts or frequencies.   
#' @usage var_plot_1000g(vcf)
#' var_plot_1000g(vcf, fraction = TRUE)
#' @param vcf plain_vcf data.frame; clinvaR-processed VCF from 1000 Genomes. 
#' @param fraction logical; if TRUE, plots fraction with a finding. If FALSE, plots total counts. 
#' @export

plot.annotated_vcf <- function(vcf, fraction, sd) {
  if (missing(fraction)) 
    fraction <- TRUE
  if (missing(sd))
    sd <- TRUE
  if (!is.logical(fraction)) {
    print("Fraction must be logical: set as TRUE")
    fraction <- TRUE
  }
  return(var_plot_1000g(vcf, fraction, sd))
}
