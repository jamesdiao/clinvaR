#' Merge ClinVar with ExAC
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage merge_clinvar_exac(clinvar, ACMG.exac)
#' @param clinvar clinvaR-processed vcf containing ClinVar data. 
#' Defaults to get_test_clinvar().  
#' @param ACMG.exac clinvaR-processed vcf containing ExAC sequencing data. 
#' Defaults to 'extdata/Supplementary_Files/ACMG_EXAC.rds'. 
#' @examples merge_clinvar_exac()
#' merge_clinvar_exac(clinvar, ACMG.exac)
#' #@export

merge_clinvar_exac <- function(clinvar, ACMG.exac) {
  if (missing(clinvar)) {
    clinvar <- get_clinvar()
  }
  if (missing(ACMG.exac)) {
    system.file("extdata", "Supplementary_Files/ACMG_EXAC.rds", package = "clinvaR") %>% 
      readRDS -> ACMG.exac
  }
  inter <- intersect(clinvar$VAR_ID[clinvar$INTERP], ACMG.exac$VAR_ID)
  merged_exac <- cbind(clinvar[(clinvar$VAR_ID %in% inter),] %>% arrange(VAR_ID), 
                       ACMG.exac %>% select(VAR_ID, contains("AF_"), GENE) %>% 
                         filter(VAR_ID %in% inter) %>% arrange(VAR_ID) %>% select(-VAR_ID)
  ) %>% select(VAR_ID, GENE, AF_EXAC, contains("AF_"), everything())
  return(merged_exac)
}



#' Merge ClinVar with gnomAD
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage merge_clinvar_gnomad(clinvar, ACMG.gnomad)
#' @param clinvar clinvaR-processed vcf containing ClinVar data. 
#' Defaults to get_test_clinvar().  
#' @param ACMG.gnomAD clinvaR-processed vcf containing gnomAD sequencing data. 
#' Defaults to 'extdata/Supplementary_Files/ACMG_GNOMAD.rds'. 
#' @examples merge_clinvar_gnomad()
#' merge_clinvar_gnomad(clinvar, ACMG.gnomad)
#' #@export

merge_clinvar_gnomad <- function(clinvar, ACMG.gnomad) {
  if (missing(clinvar)) {
    clinvar <- get_clinvar()
  }
  if (missing(ACMG.gnomad)) {
    system.file("extdata", "Supplementary_Files/ACMG_GNOMAD.rds", package = "clinvaR") %>% 
      readRDS -> ACMG.gnomad
  }
  inter <- intersect(clinvar$VAR_ID[clinvar$INTERP], ACMG.gnomad$VAR_ID)
  merged_gnomad <- cbind(clinvar[(clinvar$VAR_ID %in% inter),] %>% arrange(VAR_ID), 
                         ACMG.gnomad %>% select(VAR_ID, contains("AF_"), GENE) %>% 
                           filter(VAR_ID %in% inter) %>% arrange(VAR_ID) %>% select(-VAR_ID)
  ) %>% select(VAR_ID, GENE, AF_GNOMAD, contains("AF_"), everything())
  return(merged_gnomad)
}
