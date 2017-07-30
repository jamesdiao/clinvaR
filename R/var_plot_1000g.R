#' Compute and Plot Aggregate Allele Frequencies for 1000 Genomes
#'
#' This function takes a merged ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of allele frequencies
#' @usage var_plot_1000g(vcf)
#' var_plot_1000g(vcf, fraction = FALSE)
#' @param vcf data.frame; clinvaR-processed VCF from 1000 Genomes. Defaults to importing from genes (if present).
#' @param fraction logical; if TRUE, plots fraction with a finding. If FALSE, plots total counts. 
#' @export

var_plot_1000g <- function(vcf, fraction, sd) {
  map <- system.file("extdata", "Supplementary_Files/phase3map.txt", package = "clinvaR") %>% 
    read.table(stringsAsFactors = F, header = T) %>% as.data.frame
  pop.levels <- c("ACB", "ASW", "ESN", "GWD", "LWK", "MSL", "YRI", "CLM", "MXL", 
                  "PEL", "PUR", "CDX", "CHB", "CHS", "JPT", "KHV", "CEU", "FIN", 
                  "GBR", "IBS", "TSI", "BEB", "GIH", "ITU", "PJL", "STU")

  front_cols <- 1:(grep("HG00096",colnames(vcf))-1)
  sapply(pop.levels, function(pop) {
    keep <- c(front_cols,map$pop)==pop
    if(fraction) {
      temp <- colSums(vcf[,keep]) > 0 #average fraction with at least 1
    } else {
      temp <- colSums(vcf[,keep] > 0) #average number per person
    }
    return(c(mean(temp), sd(temp)/length(temp)))
  }) %>% t %>% as.data.frame() -> values #Number of non-reference sites across the different populations
  #colnames(values) <- c("Mean","SD")
  #values$Population <- factor(pop.levels, levels = pop.levels)
  #values$Superpopulation <- factor(super[pop.levels], levels = super.levels)
  title <- sprintf("%s in 1000 Genomes", ifelse(fraction,"Fraction","Mean Number"))
  prettyprint(values, title = title, sd = sd, ylimit = NULL, 
              xlabel = "Population", 
              ylabel = ifelse(fraction, "Fraction with at least 1 non-reference site", 
                              "Mean number of non-reference sites")
  ) -> plot.pop
  plot(plot.pop)
  return(plot.pop)
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
    theme_minimal() + theme(axis.text.x = element_text(angle = -45, hjust = 0.4)) + 
    scale_fill_discrete(name = "Superpopulation")
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
  return(plot.pop)
}


