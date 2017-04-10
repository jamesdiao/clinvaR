#' Barplot of ancestry-stratified 1000 Genomes allele frequencies
#'
#' Barplot of ancestry-stratified values 
#' @usage prettyprint(values, title = NULL, sd = TRUE, ylimit = NULL, xlabel = "Population", ylabel = NULL)
#' @examples prettyprint(values, title = title, sd = F, ylimit = NULL, 
#' xlabel = "Population", ylabel = "Fraction with at least 1 non-reference site")
#' 
prettyprint <- function(values, sd, title, xlabel, ylabel, ylimit) {
  if (missing(sd)) sd <- TRUE
  if (missing(title)) title <- NULL
  if (missing(xlabel)) xlabel <- "Population"
  if (missing(ylabel)) ylabel <- NULL
  if (missing(ylimit)) ylimit <- NULL
  colnames(values) <- c("Mean","SD")
  values$Population <- factor(pop.levels, levels = pop.levels)
  values$Superpopulation <- factor(super[pop.levels], levels = super.levels)
  
  plot.pop <- ggplot(values, aes(x=Population, y=Mean, fill = Superpopulation)) +
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