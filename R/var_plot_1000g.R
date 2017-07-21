#' Compute and Plot Aggregate Allele Frequencies for 1000 Genomes
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_1000g(pathogenic, frac)
#' @param genes character; genes of interest. Only necessary if vcf is missing. 
#' @param clinvar data.frame; clinvaR-processed VCF from clinvaR. Defaults to get_clinvar().
#' @param vcf data.frame; clinvaR-processed VCF from 1000 Genomes. Defaults to importing from genes (if present).
#' @param pathogenic logical; whether to limit frequency counts to pathogenic variants only. 
#' @param frac logical; if TRUE, plots fraction with a finding. If FALSE, plots total counts. 
#' @examples var_plot_1000g(pathogenic = T, frac = T)
#' @export

var_plot_1000g <- function(genes, clinvar, vcf, pathogenic, frac) {
  if (missing(vcf) & missing(genes)) {
    stop('Input either VCF from import_file_1000g() or a gene list.')
  }
  if (missing(vcf)) {
    vcf <- import_file_1000g(genes)
  }
  if (missing(pathogenic)) pathogenic <- TRUE
  if (missing(frac)) frac <- TRUE
  if (pathogenic){
    if (missing(clinvar)) {
      clinvar <- get_clinvar()
    }
    use_vcf <- merge_clinvar_1000g(clinvar, vcf)
  } else {
    use_vcf <- vcf
  }

  map <- system.file("extdata", "Supplementary_Files/phase3map.txt", package = "clinvaR") %>% 
    read.table(stringsAsFactors = F, header = T) %>% as.data.frame
  pop.levels <- c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", 
                  "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", 
                  "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU")

  front_cols <- 1:(grep("HG00096",colnames(use_vcf))-1)
  #recessive <- use_vcf$GENE %in% c("MUTYH","ATP7B")
  sapply(pop.levels, function(pop) {
    keep <- c(front_cols,map$pop)==pop
    temp <- ifelse(frac, 
                   colSums(use_vcf[,keep]) > 0,
                   colSums(use_vcf[,keep] > 0)
                   )
    c(mean(temp), sd(temp))
  }) %>% t %>% as.data.frame() -> values #Number of non-reference sites across the different populations
  #colnames(values) <- c("Mean","SD")
  #values$Population <- factor(pop.levels, levels = pop.levels)
  #values$Superpopulation <- factor(super[pop.levels], levels = super.levels)
  title <- sprintf("%s%s in 1000 Genomes", 
                   ifelse(pathogenic, "Pathogenic ",""), ifelse(frac,"Fraction","Mean Number"))
  prettyprint(values, title = title, sd = F, ylimit = NULL, 
              xlabel = "Population", 
              ylabel = ifelse(frac, "Fraction with at least 1 non-reference site", 
                              "Mean number of non-reference sites")
  )
}


#' Compute and Plot Aggregate Allele Frequencies for ExAC
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_exac(clinvar, ACMG.exac, pathogenic, frac)
#' @param pathogenic logical; whether to limit frequency counts to pathogenic variants only. 
#' @param frac logical; if TRUE, plots fraction with a finding. If FALSE, plots total counts. 
#' @examples var_plot_exac(clinvar, ACMG.exac, pathogenic = T, frac = T)
#' @export

var_plot_exac <- function(clinvar, ACMG.exac, pathogenic, frac) {
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
  if (frac) {
    exac_values <- data.frame(1-apply(1-exac_final, 2, prod, na.rm = T), super.levels)
  } else {
    exac_values <- data.frame(exac_final %>% colSums(na.rm = T), super.levels)
  }
  colnames(exac_values) = c("values","Superpopulation")
  ggplot(exac_values, aes(x = Superpopulation, y=values, fill = Superpopulation)) + 
    geom_bar(stat = "identity") + theme_minimal() + 
    ggtitle(sprintf("ACMG-59%s: %s in %s", 
                    ifelse(pathogenic, " Pathogenic",""), ifelse(frac,"Fraction","Mean"), dataset)) + 
    xlab("Population") + ylab("Mean No. of Non-Reference Sites") + 
    ylim(0,1.1*max(exac_values$values))
}


#' Compute and Plot Aggregate Allele Frequencies for gnomAD
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_gnomad(clinvar, ACMG.gnomad, pathogenic, frac)
#' @examples var_plot_gnomad(clinvar, ACMG.gnomad, pathogenic = T, frac = T)
#' #@export

var_plot_gnomad <- function(clinvar, ACMG.gnomad, pathogenic, frac) {
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
  if (frac) {
    exac_values <- data.frame(1-apply(1-exac_final, 2, prod, na.rm = T), super.levels)
  } else {
    exac_values <- data.frame(exac_final %>% colSums(na.rm = T), super.levels)
  }
  colnames(exac_values) = c("values","Superpopulation")
  ggplot(exac_values, aes(x = Superpopulation, y=values, fill = Superpopulation)) + 
    geom_bar(stat = "identity") + theme_minimal() + 
    ggtitle(sprintf("ACMG-59%s: %s in %s", 
                    ifelse(pathogenic, " Pathogenic",""), ifelse(frac,"Fraction","Mean"), dataset)) + 
    xlab("Population") + ylab("Mean No. of Non-Reference Sites") + 
    ylim(0,1.1*max(exac_values$values))
}


#' Formatted ggplot of ancestry-stratified 1000 Genomes allele frequencies
#'
#' Barplot of ancestry-stratified values 
#' 
#' @usage prettyprint(values, sd, title, ylimit, xlabel, ylabel)
#' ## Default method
#' prettyprint(values, title = NULL, sd = TRUE, ylimit = NULL, xlabel = "Population", ylabel = NULL)
#' @param values data frame (26 x 2); values (mean and sd) for each 1000 genomes population
#' @param sd logical; show standard deviation bars 
#' @param title character; plot title
#' @param ylimit 2 element vector; specifies bounds of y-axis
#' @param xlabel character; label for x-axis
#' @param ylabel character; label for y-axis
#' @examples 
#' prettyprint(values, title = title, sd = F, ylimit = NULL, 
#' xlabel = "Population", ylabel = "Fraction with at least 1 non-reference site")
#' @export
#' 

prettyprint <- function(values, sd, title, xlabel, ylabel, ylimit) {
  pop.levels <- c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", 
                  "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", 
                  "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU")
  super <- c("AFR", "AFR", "AFR", "AFR", "AFR", "AFR", "AFR", "AMR", "AMR", 
             "AMR", "AMR", "EAS", "EAS", "EAS", "EAS", "EAS", "EUR", "EUR", 
             "EUR", "EUR", "EUR", "SAS", "SAS", "SAS", "SAS", "SAS") %>% setNames(pop.levels)
  super.levels <- unique(super)
  
  if (missing(sd)) sd <- TRUE
  if (missing(title)) title <- NULL
  if (missing(xlabel)) xlabel <- "Population"
  if (missing(ylabel)) ylabel <- NULL
  if (missing(ylimit)) ylimit <- NULL
  colnames(values) <- c("Mean","SD")
  values$Population <- factor(pop.levels, levels = pop.levels)
  values$Superpopulation <- factor(super[pop.levels], levels = super.levels)
  
  plot.pop <- ggplot(values, aes(x=values$Population, y=Mean, fill = values$Superpopulation)) +
    geom_bar(stat = "identity") + ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
    theme_minimal() + theme(axis.text.x = element_text(angle = -45, hjust = 0.4))
  if (sd) {
    if (min(values$Mean - values$SD)<0)
      plot.pop <- plot.pop + geom_errorbar(aes(
        ymin = pmax(0,values$Mean - values$SD), 
        ymax = values$Mean + values$SD, width = 0.5))
    else 
      plot.pop <- plot.pop + geom_errorbar(aes(ymin = values$Mean - values$SD, 
                                               ymax = values$Mean + values$SD, width = 0.5))
  } else {values$SD = 0}
  if (length(ylimit)==2)
    plot.pop <- plot.pop + ylim(ylimit[1],ylimit[2])
  else
    plot.pop <- plot.pop + ylim(0, 1.1*max(values$Mean + values$SD))
  plot.pop
}


