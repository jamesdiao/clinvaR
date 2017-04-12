#' Compute Aggregate Allele Frequencies by ACMG Condition
#'
#' @usage getAlleleFreq(input, ind, dataset, method)
#' @examples freq_1000g.count.gene <- getAlleleFreq(input = merged_1000g, ind = F, dataset = "1000G", method = "Gene")
#' freq_1000g.calc.gene <- getAlleleFreq(input = merged_1000g, ind = T, dataset = "1000G", method = "Gene")
#' freq_gnomad.calc.gene <- getAlleleFreq(input = merged_gnomad, ind = T, dataset = "GNOMAD", method = "Gene")
#' freq_exac.calc.gene <- getAlleleFreq(input = merged_exac, ind = T, dataset = "EXAC",method = "Gene")
#' @export

getAlleleFreq <- function(input, ind, dataset, method) {
  #if (!("CLNSIG" %in% colnames(input)))
  #  input$CLNSIG <- rep(5, nrow(input))
  if (toupper(method) == "GENE") {
    expand = F; search_in = "GENE"; search_col = "Gene"; 
    report.use = report.gene; inherit.use = inheritance.gene
  }
  if (toupper(method) == "MIM") {
    expand = F; search_in = "CLNDSDBID"; search_col = "MIM"; 
    report.use = report.mim; inherit.use = inheritance.mim
  }
  if (toupper(method) == "MEDGEN") {
    expand = F; search_in = "CLNDSDBID"; search_col = "MedGen"; 
    report.use = report.medgen; inherit.use = inheritance.medgen
  }
  if (expand)
    input <- input %>% filter(list_len(input$CLNSIG) == list_len(input$CLNDSDBID)) %>% 
      unnest %>% select(VAR_ID, GENE, contains("AF"), CHROM, POS, ID, REF, ALT, QUAL, 
                        FILTER, INTERP, contains("CLN"), INFO, everything())
  KP_only <- grepl(5,input$CLNSIG)
  data.frame(ACMG_Lit_Full$Short_Name, 
             sapply(1:nrow(ACMG_Lit_Full), function(index) {
               item.vec.split <- expand_pipes(ACMG_Lit_Full[[search_col]][index]) %>% unique
               sapply(item.vec.split, function(item) {
                 loc <- grepl(item,input[,search_in], ignore.case = T)
                 if (!grepl("EP",report.use[[item]])) #If we're only taking KP
                   loc <- loc & KP_only #Take all relevant genes
                 sapply(c(dataset,super.levels), function(superpop) {
                   if (ind) {
                     return(aggregateCalc(input, superpop, item, dataset, loc, inherit.use))
                   } else {
                     return(aggregateCount(input, superpop, item, dataset, loc, inherit.use))
                   }
                 })
               }) %>% apply(1, function(row) 1-prod(1-row)) %>% 
                 setNames(sprintf("AF_%s%s",toupper(dataset), c("",paste0("_",super.levels))))
             }) %>% t %>% tbl_df
  )
}