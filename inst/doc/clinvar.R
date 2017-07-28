## ----setup, include = F--------------------------------------------------
knitr::opts_knit$set(root.dir = ".");
knitr::opts_chunk$set(echo = T, eval = T, cache = F, warning = F, message = F)

## ---- echo = F-----------------------------------------------------------
pkg_list <- c("remotes","pander","ggplot2","tibble","tidyr","dplyr", "stringr")
installed <- installed.packages()[,"Package"]
is.installed <- pkg_list %in% installed
if (!all(is.installed))
  install.packages(pkg_list[!is.installed])
print(names(sapply(pkg_list, require, character.only = T)), quote = F)
if (!("clinvaR" %in% installed))
  remotes::install_github('jamesdiao/clinvaR')
require(clinvaR)

## ------------------------------------------------------------------------
gene_panel <- c("ACTC1", "ACTN2", "CSRP3", "GLA", "LAMP2", "MYBPC3", "MYH7", 
"MYL2", "MYL3", "MYOZ2", "NEXN", "PLN", "PRKAG2", "PTPN11", "RAF1", 
"TNNC1", "TNNI3", "TNNT2", "TPM1", "TTR")

## ------------------------------------------------------------------------
path <- system.file("extdata/MacArthur_Gene_Lists/lists/lmm_hcm.tsv", package = "clinvaR")
print(path)

## ------------------------------------------------------------------------
gene_panel <- get_genes(path)
print(gene_panel)

## ------------------------------------------------------------------------
gene_panel <- get_genes('lmm_hcm.tsv')
print(gene_panel)

## ------------------------------------------------------------------------
#omim_table <- read.table('clinvaR/storage/gene_lists/other_data/omim.full.tsv', 
#                         sep = '\t', fill = T, header = T)
#search_terms <- c("Cardiomyopathy")

## ------------------------------------------------------------------------
#download_output <- download_1000g(genes = gene_panel)
#final_vcf <- import_file_1000g(genes = gene_panel)

## ------------------------------------------------------------------------
get_date_list()

## ------------------------------------------------------------------------
clinvar <- get_clinvar("2017-07-21")
clinvar[1:4,] #%>% select(-VAR_ID, -CLNDBN)

## ------------------------------------------------------------------------
clinvar <- get_clinvar()
clinvar[1:4,] #%>% select(-VAR_ID, -CLNDBN)

## ---- cache = T----------------------------------------------------------
#merged_1000g <- merge_clinvar_1000g(clinvar = get_clinvar(Sys.Date()), 
#                    vcf = import_file_1000g(gene_panel))
#merged_1000g[1:5,1:5]

## ------------------------------------------------------------------------
#clinvar$CHROM %>% paste0(str_pad(clinvar$POS, 10, side = 'left', pad = "0")) %>% as.numeric %>% 
#  between(start, stop) -> filter_condition
#clinvar <- filter(clinvar, filter_condition)

