#' Compute and Plot Aggregate Allele Frequencies for 1000 Genomes
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_1000g(pathogenic, frac)
#' @examples var_plot_1000g(pathogenic = T, frac = T)
#' @export

var_plot_1000g <- function(clinvar, ACMG.1000g, pathogenic, frac) {
  merged_1000g <- merge_clinvar_1000g(clinvar, ACMG.1000g)
  if (pathogenic){
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
  colnames(values) <- c("Mean","SD")
  values$Population <- factor(pop.levels, levels = pop.levels)
  values$Superpopulation <- factor(super[pop.levels], levels = super.levels)
  title <- sprintf("ACMG-59%s: %s in 1000 Genomes", 
                   ifelse(pathogenic, " Pathogenic",""), ifelse(frac,"Fraction","Mean"))
  prettyprint(values, title = title, sd = F, ylimit = NULL, 
              xlabel = "Population", 
              ylabel = ifelse(frac, "Fraction with at least 1 non-reference site", 
                              "Mean No. of Non-Reference Sites")
  )
}