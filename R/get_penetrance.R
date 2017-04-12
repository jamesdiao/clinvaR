#' Compute and Plot Penetrance Ranges by ACMG Condition
#'
#' @usage get_penetrance(ah_low, ah_high, dataset)
#' @examples pen_gnomad <- get_penetrance(ah_low = 0.01, ah_high = 1, dataset = "gnomAD")
#' pen_1000g <-  get_penetrance(ah_low = 0.01, ah_high = 1, dataset = "1000 Genomes")
#' @export

get_penetrance <- function(ah_low, ah_high, dataset) {
  # Map of disease name to disease tags
  if (dataset == "1000 Genomes")
    named.freqs <- allele.freq$COUNT_1000G %>% setNames(abbrev)
  if (dataset == "gnomAD")
    named.freqs <- allele.freq$CALC_GNOMAD %>% setNames(abbrev)
  if (dataset == "ExAC")
    named.freqs <- allele.freq$CALC_EXAC %>% setNames(abbrev)
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

