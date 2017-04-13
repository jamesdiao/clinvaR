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
  if (missing(dataset)) dataset <- 'GNOMAD'
  if (missing(input)) input <- eval(parse(text=sprintf('merge_clinvar_%s()', tolower(dataset)) ))
  if (missing(ind)) ind <- TRUE
  if (missing(method)) method <- "Gene"
  
  super.levels <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  map <- system.file("extdata", "Supplementary_Files/phase3map.txt", package = "clinvaR") %>% 
    read.table(stringsAsFactors = F, header = T) %>% as.data.frame
  sample_size <- list()
  sample_size$Cohort_1000G <- table(map$super_pop) %>% 
    as.numeric %>% setNames(super.levels) %>% c("1000G"=2504)
  sample_size$Cohort_EXAC <- c("AFR"=5203,"AMR"=5789,"EAS"=4327,
                               "EUR"=3307+33370,"SAS"=8256,"EXAC"=60252)
  sample_size$Cohort_GNOMAD <- c("AFR"=12942,"AMR"=18237,"EAS"=9472,
                                 "EUR"=5081+13046+63416,"SAS"=15450,"GNOMAD"=137644)
  
  ACMG_Lit_Full <- system.file("extdata", "Supplementary_Files/Literature_Prevalence_Estimates.csv", package = "clinvaR") %>% 
    read.csv(header = TRUE, stringsAsFactors = F, na.strings = "\\N") 
  expand_pipes <- function(item) { strsplit(item, "|", fixed = T) %>% unlist }
  gene.list <- expand_pipes(ACMG_Lit_Full$Gene)
  MIM.list <- expand_pipes(ACMG_Lit_Full$MIM)
  MG.list <- expand_pipes(ACMG_Lit_Full$MedGen)
  report <- list()
    report$all <- ACMG_Lit_Full$Variants_to_report %>% expand_pipes()
    report$gene <- setNames(report$all, gene.list)[!duplicated(gene.list)]
    report$mim <- setNames(report$all, MIM.list)[!duplicated(MIM.list)]
    report$medgen <- setNames(report$all, MG.list)[!duplicated(MG.list)]
  inheritance <- list()
    inheritance$all <- ACMG_Lit_Full$Inheritance %>% expand_pipes()
    inheritance$gene <- setNames(inheritance$all, gene.list)[!duplicated(gene.list)]
    inheritance$mim <- setNames(inheritance$all, MIM.list)[!duplicated(MIM.list)]
    inheritance$medgen <- setNames(inheritance$all, MG.list)[!duplicated(MG.list)]
  
  list_len <- function(list_in) { unlist(lapply(list_in, length)) }
  if (toupper(method) == "GENE") {
    expand = F; search_in = "GENE"; search_col = "Gene"; 
    report.use = report$gene; inherit.use = inheritance$gene
  }
  if (toupper(method) == "MIM") {
    expand = F; search_in = "CLNDSDBID"; search_col = "MIM"; 
    report.use = report$mim; inherit.use = inheritance$mim
  }
  if (toupper(method) == "MEDGEN") {
    expand = F; search_in = "CLNDSDBID"; search_col = "MedGen"; 
    report.use = report$medgen; inherit.use = inheritance$medgen
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
                     return(aggregateCalc(input, superpop, item, dataset, loc, inherit.use, sample_size))
                   } else {
                     return(aggregateCount(input, superpop, item, dataset, loc, inherit.use, sample_size))
                   }
                 })
               }) %>% apply(1, function(row) 1-prod(1-row)) %>% 
                 setNames(sprintf("AF_%s%s",toupper(dataset), c("",paste0("_",super.levels))))
             }) %>% t %>% tbl_df
  )
}