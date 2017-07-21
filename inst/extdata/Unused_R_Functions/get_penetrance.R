#' Compute and Plot Penetrance Ranges by ACMG Condition
#'
#' @usage get_penetrance(ah_low, ah_high, dataset, alleleFreq)
#' @param ah_low numeric; lower allelic heterogeneity bound.
#' @param ah_high numeric; higher allelic heterogeneity bound.
#' @param dataset character; dataset name. Choose from 1000 Genomes, ExAC, or gnomAD (not case sensitive)
#' @param alleleFreq data frame (31 x n); contains the columns 'AF_[dataset]_[population]' 
#' for 5 superpopulations and the general population. Typically an output from getAlleleFreq().
#' @examples 
#' pen_gnomad <- get_penetrance(ah_low = 0.01, ah_high = 1, dataset = "gnomAD", alleleFreq = freq_gnomad.calc.gene)
#' pen_1000g <-  get_penetrance(ah_low = 0.01, ah_high = 1, dataset = "1000 Genomes", alleleFreq = freq_1000g.count.gene)
#' #@export

get_penetrance <- function(ah_low, ah_high, dataset, alleleFreq) {
  # Map of disease name to disease tags
  ACMG_Lit_Full <- system.file("extdata", "Supplementary_Files/Literature_Prevalence_Estimates.csv", package = "clinvaR") %>% 
    read.csv(header = TRUE, stringsAsFactors = F, na.strings = "\\N") 
  ACMG_Lit <- ACMG_Lit_Full %>% filter(Evaluate)
  abbrev <- ACMG_Lit$Short_Name
  if (nrow(alleleFreq)==nrow(ACMG_Lit_Full))
    alleleFreq <- alleleFreq[ACMG_Lit_Full$Evaluate,]
  inv.prev <- ACMG_Lit$Inverse_Prevalence %>% as.numeric %>% setNames(abbrev)
  acmg_ah <- ACMG_Lit$Case_Allele_Frequency %>% as.numeric
  
  if (toupper(dataset) == "1000 GENOMES")
    named.freqs <- alleleFreq$AF_1000G %>% setNames(abbrev)
  if (toupper(dataset) == "GNOMAD")
    named.freqs <- alleleFreq$AF_GNOMAD %>% setNames(abbrev)
  if (toupper(dataset) == "EXAC")
    named.freqs <- alleleFreq$AF_EXAC %>% setNames(abbrev)
  named.prev <- 1/inv.prev %>% setNames(abbrev)
  # Repeats allow for correct quartile calculations
  #point estimate set to arithmetic mean
  allelic.het <- c(ah_low, ah_low, mean(c(ah_low, ah_high)) %>% signif(3), ah_high, ah_high) %>% 
    rep(nrow(ACMG_Lit)) %>% matrix(nrow = nrow(ACMG_Lit), byrow = T)
  allelic.het[,3] <- acmg_ah
  #Functions to transform data points with disease_af = 0
  set_to_na <- function(row) { replace(row, is.infinite(row),NA) %>% pmin(1)}
  # Matrix of penetrance values for allelic het range, capped at 1
  prev_freq <- named.prev/named.freqs %>% rep(5) %>% matrix(nrow = nrow(ACMG_Lit))
  penetrance <- apply(allelic.het * prev_freq, 1, set_to_na) %>% as.data.frame
  # Take row 5 to sort by max, take colSums to sort by overall
  # Break ties by rows 3 and 1 (mean and low)
  ord <- order(penetrance[5,], penetrance[3,], penetrance[1,], decreasing = T)
  # replicate each element n times to create labels
  penetrance_data <- data.frame("Penetrance" = penetrance %>% unlist, "Disease" =
                                  factor(sapply(abbrev, function(x) rep(x,5)) %>% as.vector,
                                         levels = abbrev[ord]) )
  colormap <- rep('black', length(abbrev))
  num_na <- mean(is.na(penetrance_data$Penetrance))*length(colormap)
  if (num_na != 0)
    colormap[1:num_na] <- 'gray60'
  colormap <- rev(colormap)
  plot(
    ggplot(aes(x=Disease, y=Penetrance), data = penetrance_data) + 
      geom_boxplot(position = 'identity', coef = 0, na.rm = T) + coord_flip() + 
      ggtitle(sprintf("%s: Barplot of Min/Point/Max Penetrance", dataset)) + 
      theme(axis.text.y = element_text(color = colormap)) + ylim(0,1)
    #annotate("text", x = which(colormap == 'red'), y = 0, 
    #         label = "No Allele Frequency Data", hjust = 0, size = 3)
  )
  penetrance_data
}

