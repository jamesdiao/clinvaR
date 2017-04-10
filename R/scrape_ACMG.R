#' Scrape ACMG Table from ClinVar
#'
#' This function scrapes the table from 'https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg' and cleans it up. 
#' @usage scrape_ACMG(scrape)
#' @keywords ACMG

scrape_ACMG <- function(scrape) {
  if (missing(scrape)) {
    scrape <- TRUE 
  }
  if (scrape) {
    ACMG.page <- scrape(url ="https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/")[[1]]
    ACMG.table <- readHTMLTable(ACMG.page, stringsAsFactors = F, header = T)[[1]]
    colnames(ACMG.table) <- c("Disease", "MedGen","Gene","Variations_Link")
    ### Formatting correction: separating merged gene blocks
    badrow <- which(apply(ACMG.table, 1, function(row) !any(grepl("ClinVar", row))))
    ACMG.table[badrow,"Gene"] <- ACMG.table[badrow-1,"Gene"]
    ### Formatting corrections: sliding
    mismatch <- 0
    while(any(ACMG.table[,"Gene"] == "ClinVar")) {
      mismatch <- which(ACMG.table[,"MedGen"]!="MedGen")
      ACMG.table[mismatch,2:3] <- ACMG.table[mismatch,1:2]
      for (row in mismatch) { ACMG.table[row,"Disease"] <- ACMG.table[row-1, "Disease"] }
    }
    ACMG.table %>% select(Disease, Gene) %>% 
      separate(col = Disease, into = c("Disease_Name","Disease_MIM"), sep = " \\(MIM ") %>%
      separate(col = Gene, into = c("Gene_Name","Gene_MIM"), sep = " \\(MIM ") %>% 
      mutate(Disease_MIM = strsplit(Disease_MIM, "\\)") %>% sapply("[",1)) %>%
      mutate(Gene_MIM = strsplit(Gene_MIM, "\\)") %>% sapply("[",1))
    
  } else {
    ACMG.table <- read.table(file = "Supplementary_Files/ACMG_SF_v2.0.txt", 
                             stringsAsFactors = F, sep = "\t", header = T) %>%
      mutate(Gene = strsplit(Gene, "|", fixed = T)) %>% 
      mutate(MIM_gene = strsplit(MIM_gene, "|", fixed = T)) %>% 
      mutate(Inheritance = strsplit(Inheritance, "|", fixed = T)) %>% 
      mutate(Variants_to_report = strsplit(Variants_to_report, "|", fixed = T)) %>% 
      unnest()
  }

}