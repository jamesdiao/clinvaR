## ----setup, include = F--------------------------------------------------
knitr::opts_knit$set(root.dir = ".");
knitr::opts_chunk$set(echo = T, eval = T, cache = T, warning = F, message = F)

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
hcm_panel <- c("ACTC1", "ACTN2", "CSRP3", "GLA", "LAMP2", "MYBPC3", "MYH7", 
"MYL2", "MYL3", "MYOZ2", "NEXN", "PLN", "PRKAG2", "PTPN11", "RAF1", 
"TNNC1", "TNNI3", "TNNT2", "TPM1", "TTR")

## ------------------------------------------------------------------------
path <- system.file("extdata/MacArthur_Gene_Lists/lists/lmm_hcm.tsv", package = "clinvaR")
print(path)
hcm_panel <- get_genes(path)
print(hcm_panel)

## ------------------------------------------------------------------------
get_genes()
hcm_panel <- get_genes('lmm_hcm.tsv')
print(hcm_panel)

## ---- echo = F, include = F----------------------------------------------
#omim_table <- read.table('clinvaR/storage/gene_lists/other_data/omim.full.tsv', 
#                         sep = '\t', fill = T, header = T)
#search_terms <- c("Cardiomyopathy")

## ------------------------------------------------------------------------
#hcm_download_log <- download_1000g(genes = hcm_panel)
vcf <- import_file_1000g(genes = hcm_panel)
vcf[1:3,1:8]

## ------------------------------------------------------------------------
get_date_list()

## ------------------------------------------------------------------------
newest_clinvar <- get_clinvar(Sys.Date())
newest_clinvar[1:4,] #%>% select(-VAR_ID, -CLNDBN)

## ---- cache = T----------------------------------------------------------
hcm_vcf <- annotate_1000g(vcf = vcf, clinvar = newest_clinvar, conflicts = TRUE)
hcm_vcf[1:3,1:12]

## ---- fig.width = 6.5, fig.height = 4.5----------------------------------
plot(hcm_vcf, fraction = TRUE)
plot(hcm_vcf, fraction = FALSE)

## ---- fig.width = 6.5, fig.height = 4.5----------------------------------
plot(vcf, fraction = FALSE)

