#' Merge ClinVar with 1000 Genomes
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage merge_clinvar_1000g(clinvar, ACMG.1000g)
#' @examples merge_clinvar_1000g(clinvar, ACMG.1000g)
#' @export

merge_clinvar_1000g <- function(clinvar, ACMG.1000g) {
  if (missing(ACMG.1000g)) {
    system.file("extdata", "Supplementary_Files/ACMG_1000G.rds", package = "clinvaR") %>% 
      readRDS -> ACMG.1000g
  }
  inter <- intersect(clinvar$VAR_ID[clinvar$INTERP], ACMG.1000g$VAR_ID)
  clinvar_merged <- clinvar[(clinvar$VAR_ID %in% inter),] %>% arrange(VAR_ID)
  ACMG_merged <- ACMG.1000g[ACMG.1000g$VAR_ID %in% inter,] %>% arrange(VAR_ID)
  front_cols <- 1:(grep("HG00096",colnames(ACMG.1000g))-1)
  super.levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  ACMG_merged <- cbind(ACMG_merged[,c("GENE",sprintf("AF_1000G%s", c("",paste0("_",super.levels))))], 
        clinvar_merged,ACMG_merged[,-front_cols])
  return(ACMG_merged)
}

#' Merge ClinVar with ExAC
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage merge_clinvar_exac(clinvar, ACMG.1000g)
#' @examples merge_clinvar_exac(clinvar, ACMG.1000g)
#' @export
# ADD ACMG_EXAC.rds
merge_clinvar_exac <- function(clinvar, ACMG.exac) {
  if (missing(ACMG.exac)) {
    system.file("extdata", "Supplementary_Files/ACMG_1000G.rds", package = "clinvaR") %>% 
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
#' @usage merge_clinvar_gnomad(clinvar, ACMG.1000g)
#' @examples merge_clinvar_gnomad(clinvar, ACMG.1000g)
#' @export

merge_clinvar_gnomad <- function() {
  inter <- intersect(clinvar$VAR_ID[clinvar$INTERP], ACMG.gnomad$VAR_ID)
  merged_gnomad <- cbind(clinvar[(clinvar$VAR_ID %in% inter),] %>% arrange(VAR_ID), 
                         ACMG.gnomad %>% select(VAR_ID, contains("AF_"), GENE) %>% 
                           filter(VAR_ID %in% inter) %>% arrange(VAR_ID) %>% select(-VAR_ID)
  ) %>% select(VAR_ID, GENE, AF_GNOMAD, contains("AF_"), everything())
  return(merged_gnomad)
}

