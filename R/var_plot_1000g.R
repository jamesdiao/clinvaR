#' Compute and Plot Aggregate Allele Frequencies for 1000 Genomes
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_1000g(pathogenic, frac)
#' @examples var_plot_1000g(pathogenic = T, frac = T)
#' @export

var_plot_1000g <- function(clinvar, ACMG.1000g, pathogenic, frac) {
  if (missing(clinvar)) {
    clinvar <- get_test_clinvar()
  }
  if (missing(ACMG.1000g)) {
    system.file("extdata", "Supplementary_Files/ACMG_1000G.rds", package = "clinvaR") %>% 
      readRDS -> ACMG.1000g
  }
  map <- system.file("extdata", "Supplementary_Files/phase3map.txt", package = "clinvaR") %>% 
    read.table(stringsAsFactors = F, header = T) %>% as.data.frame
  pop.levels <- c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", 
                  "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", 
                  "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU")
  if (pathogenic){
    merged_1000g <- merge_clinvar_1000g(clinvar, ACMG.1000g)
    KP <- sapply(merged_1000g$CLNSIG, function(x) 5 %in% x)
    KP_only <- c("RET","PRKAG2","MYH7","TNNI3","TPM1","MYL3","CACNA1S",
                 "DSP","MYL2","APOB","PCSK9","RYR1","RYR2","SDHAF2","ACTC1")
    ACMG.data <- merged_1000g[!((merged_1000g$GENE %in% KP_only) & (!KP)),]
  } else {
    ACMG.data <- ACMG.1000g
  }
  front_cols <- 1:(grep("HG00096",colnames(ACMG.data))-1)
  recessive <- ACMG.data$GENE %in% c("MUTYH","ATP7B")
  sapply(pop.levels, function(pop) {
    keep <- c(front_cols,map$pop)==pop
    if (frac) {
      temp <- (colSums(ACMG.data[!recessive,keep])>0) + 
        (colSums(ACMG.data[ recessive,keep])>1)
    } else {
      #Counts the number of non-reference sites in a gene
      temp <- colSums(ACMG.data[,keep]>0)
    }
    c(mean(temp), sd(temp))
  }) %>% t %>% tbl_df -> values #Number of non-reference sites across the different populations
  #colnames(values) <- c("Mean","SD")
  #values$Population <- factor(pop.levels, levels = pop.levels)
  #values$Superpopulation <- factor(super[pop.levels], levels = super.levels)
  title <- sprintf("ACMG-59%s: %s in 1000 Genomes", 
                   ifelse(pathogenic, " Pathogenic",""), ifelse(frac,"Fraction","Mean"))
  prettyprint(values, title = title, sd = F, ylimit = NULL, 
              xlabel = "Population", 
              ylabel = ifelse(frac, "Fraction with at least 1 non-reference site", 
                              "Mean No. of Non-Reference Sites")
  )
}


#' Compute and Plot Aggregate Allele Frequencies for ExAC
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_exac(clinvar, ACMG.exac, pathogenic, frac)
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
#' @export

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