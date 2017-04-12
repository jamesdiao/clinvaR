#' Compute and Plot Ancestry-Stratified Max/Median Penetrance by ACMG Condition
#'
#' @usage ancestry_penetrance(ah_low, ah_high, dataset, range, position, af_type)
#' @examples ancestry_pen <- ancestry_penetrance(ah_low = 0.01, ah_high = 1, dataset = "gnomAD", 
#' range = 5, position = "Max", af_type = "gene")
#' @export

ancestry_penetrance <- function(ah_low, ah_high, dataset, range, position, af_type) {
  pos <- replace(c(F,F,F,F,F), ifelse(position == "Max", 5, 3), T)
  sapply(c("Total",super.levels), function(superpop){
    # Map of disease name to disease tags
    find <- paste0("AF_", toupper(dataset))
    named.freqs <- eval(parse(text=sprintf("freq_%s.calc.%s", tolower(dataset), af_type)))[ACMG_Lit_Full$Evaluate,]
    if (superpop != "Total") 
      find <- paste(find, superpop, sep = "_")
    named.freqs <- named.freqs[,find] %>% unlist %>% setNames(abbrev)
    allelic.het <- c(ah_low, ah_low, mean(c(ah_low, ah_high)) %>% signif(3), ah_high, ah_high) %>% 
      rep(nrow(ACMG_Lit)) %>% matrix(nrow = nrow(ACMG_Lit), byrow = T)
    allelic.het[,3] <- acmg_ah
    #make_prev_ranges(range) -> prev_final
    # Matrix of penetrance values for allelic het range, capped at 1
    set_to_na <- function(row) { replace(row, is.infinite(row),NA) %>% pmin(1)}
    apply(allelic.het / ACMG_Lit$Inverse_Prevalence / named.freqs, 1, set_to_na) %>% unlist
  }) %>% as.data.frame -> penetrance_init
  ### Quantifying ancestral variation
  #temp <- penetrance_init[pos,]
  #sort(apply(temp, 1, max, na.rm = T) - apply(temp, 1, min, na.rm = T))
  #sort(apply(temp, 1, max, na.rm = T) / apply(temp, 1, min, na.rm = T))
  # Take column 5 to sort by max, take rowSums to sort by overall
  # Break ties by rows 3 and 1 (mean and low)
  mat_pen <- matrix(penetrance_init$Total, ncol = 5, byrow = T)
  ord <- order(mat_pen[,5], mat_pen[,3], mat_pen[,1], decreasing = T)
  output_data <- data.frame(penetrance_init, 
                            "Disease" = factor(sapply(abbrev, 
                                                      function(x) rep(x,5)) %>% as.vector,
                                               levels = abbrev[ord]) ) 
  #m <- list(l = 170, r = 220, b = -50, t = 50, pad = 5)
  #  heatmap <- plot_ly(
  #    x = factor(c(super.levels,"Total"), levels = c("Total",super.levels)),
  #    y = factor(sapply(abbrev, function(x) rep(x,5)) %>% as.vector, levels = abbrev[ord]),
  #    z = penetrance_init[pos,][ord,] %>% as.matrix, type = "heatmap", height = 700
  #  ) %>% layout(autosize = T, margin = m); heatmap
  # Star/Radar Plot
  temp <- output_data[pos,] %>% select(-Disease, -Total)
  rownames(temp) <- abbrev
  order(apply(temp, 1, function(row) var(row, na.rm = T)), decreasing = T)[1:10] -> wanted
  col <- 3
  stars(temp[wanted,], scale = F, full = F, len = 1, nrow = ceiling(10/col), ncol = col, 
        flip.labels = F, key.loc = c(2*col,2), 
        main = sprintf("Radar Plot: %s Penetrance by Ancestry (%s)", position, dataset),
        draw.segments = T, col.segments = c('red','yellow','green','blue','purple'))
  print("These are the top 10 diseases by summed allele frequencies. NULL values are not plotted.", quote = F)
  print("Each radius is proportional to the penetrance of the disease in the given population.", quote = F)
  # Barplot
  penetrance_data <- gather(output_data, Subset, Penetrance, -Disease)
  plot(ggplot(aes(x=Disease, y=Penetrance), data = penetrance_data) + 
         geom_boxplot(position = 'identity', coef = 0, na.rm = T) + coord_flip() + 
         facet_wrap(~Subset, ncol = 2) + ggtitle(sprintf("Barplot: Penetrance by Ancestry (%s)", dataset)) + 
         theme(axis.text.y=element_text(size=6), axis.text.x = element_text(angle = -20, hjust = 0.4))
  )
  # Heatmap
  plot(ggplot(aes(x=Disease, y = Subset), data = penetrance_data[pos,]) + coord_flip() + 
         geom_tile(aes(fill = Penetrance), color = 'white') + xlab("Disease") + ylab("Ancestry") +
         scale_fill_gradient(low='white',high = 'darkblue', na.value = "grey50",
                             breaks=c(0,0.25,0.5,0.75,1), labels=c("0","0.25","0.50","0.75","1.00"), limits =c(0,1)) + 
         ggtitle(sprintf("Heatmap: %s Penetrance by Ancestry (%s)", position, dataset)) + 
         theme_minimal() + theme(axis.ticks = element_blank()) + 
         annotate("segment", y=c(0.5,5.5,6.5), yend=c(0.5,5.5,6.5), 
                  x=0.5, xend = length(abbrev)+0.5) +
         annotate("segment", y=0.5, yend=6.5, 
                  x=c(0.5,length(abbrev)+0.5), 
                  xend = c(0.5,length(abbrev)+0.5))
  )
  output_data[pos,]
}