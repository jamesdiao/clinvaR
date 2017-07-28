#' Compute and Plot Aggregate Allele Frequencies for ExAC
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_exac(clinvar, ACMG.exac, pathogenic, fraction)
#' @param pathogenic logical; whether to limit frequency counts to pathogenic variants only. 
#' @param fraction logical; if TRUE, plots fraction with a finding. If FALSE, plots total counts. 
#' @examples var_plot_exac(clinvar, ACMG.exac, pathogenic = T, fraction = T)
#' @export

var_plot_exac <- function(clinvar, ACMG.exac, pathogenic, fraction) {
  dataset <- 'ExAC'
  super.levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  if (missing(clinvar)) {
    clinvar <- get_test_clinvar()
  }
  if (missing(ACMG.exac)) {
    system.file("extdata", "Supplementary_Files/ACMG_EXAC.rds", package = "clinvaR") %>% 
      readRDS -> ACMG.exac
  }
  if (pathogenic) {
    ACMG.data <- merge_clinvar_exac(clinvar, ACMG.exac)
    KP <- sapply(ACMG.data$CLNSIG, function(x) 5 %in% x) 
    KP_only <- c("RET","PRKAG2","MYH7","TNNI3","TPM1","MYL3","CACNA1S",
                 "DSP","MYL2","APOB","PCSK9","RYR1","RYR2","SDHAF2","ACTC1")
    ACMG.data <- ACMG.data[!((ACMG.data$GENE %in% KP_only) & (!KP)),]
  } else {
    ACMG.data <- ACMG.exac
  }
  recessive <- ACMG.data$GENE %in% c("MUTYH","ATP7B")
  exac_prob <- ACMG.data[,sprintf("AF_%s_%s", toupper(dataset), super.levels)]
  exac_final <- 1-(1-exac_prob)^2
  exac_final[recessive,] <- (exac_prob[recessive,])^2
  if (fraction) {
    exac_values <- data.frame(1-apply(1-exac_final, 2, prod, na.rm = T), super.levels)
  } else {
    exac_values <- data.frame(exac_final %>% colSums(na.rm = T), super.levels)
  }
  colnames(exac_values) = c("values","Superpopulation")
  ggplot(exac_values, aes(x = Superpopulation, y=values, fill = Superpopulation)) + 
    geom_bar(stat = "identity") + theme_minimal() + 
    ggtitle(sprintf("ACMG-59%s: %s in %s", 
                    ifelse(pathogenic, " Pathogenic",""), ifelse(fraction,"Fraction","Mean"), dataset)) + 
    xlab("Population") + ylab("Mean No. of Non-Reference Sites") + 
    ylim(0,1.1*max(exac_values$values))
}


#' Compute and Plot Aggregate Allele Frequencies for gnomAD
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_gnomad(clinvar, ACMG.gnomad, pathogenic, fraction)
#' @examples var_plot_gnomad(clinvar, ACMG.gnomad, pathogenic = T, fraction = T)
#' #@export

var_plot_gnomad <- function(clinvar, ACMG.gnomad, pathogenic, fraction) {
  dataset <- 'gnomAD'
  super.levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  if (missing(clinvar)) {
    clinvar <- get_test_clinvar()
  }
  if (missing(ACMG.gnomad)) {
    system.file("extdata", "Supplementary_Files/ACMG_GNOMAD.rds", package = "clinvaR") %>% 
      readRDS -> ACMG.gnomad
  }
  if (pathogenic) {
    ACMG.data <- merge_clinvar_gnomad(clinvar, ACMG.gnomad)
    KP <- sapply(ACMG.data$CLNSIG, function(x) 5 %in% x) 
    KP_only <- c("RET","PRKAG2","MYH7","TNNI3","TPM1","MYL3","CACNA1S",
                 "DSP","MYL2","APOB","PCSK9","RYR1","RYR2","SDHAF2","ACTC1")
    ACMG.data <- ACMG.data[!((ACMG.data$GENE %in% KP_only) & (!KP)),]
  } else {
    ACMG.data <- ACMG.gnomad
  }
  recessive <- ACMG.data$GENE %in% c("MUTYH","ATP7B")
  exac_prob <- ACMG.data[,sprintf("AF_%s_%s", toupper(dataset), super.levels)]
  exac_final <- 1-(1-exac_prob)^2
  exac_final[recessive,] <- (exac_prob[recessive,])^2
  if (fraction) {
    exac_values <- data.frame(1-apply(1-exac_final, 2, prod, na.rm = T), super.levels)
  } else {
    exac_values <- data.frame(exac_final %>% colSums(na.rm = T), super.levels)
  }
  colnames(exac_values) = c("values","Superpopulation")
  ggplot(exac_values, aes(x = Superpopulation, y=values, fill = Superpopulation)) + 
    geom_bar(stat = "identity") + theme_minimal() + 
    ggtitle(sprintf("ACMG-59%s: %s in %s", 
                    ifelse(pathogenic, " Pathogenic",""), ifelse(fraction,"Fraction","Mean"), dataset)) + 
    xlab("Population") + ylab("Mean No. of Non-Reference Sites") + 
    ylim(0,1.1*max(exac_values$values))
}