#' Imports downloaded VCF from ExAC or gnomAD
#'
#' This function imports a gene's VCF and extracts attributes into a 
#' specified format for other clinvaR analyses
#' 
#' @usage import_file_exac(gene, dataset)
#' @examples import_file_exac('BRCA2', 'gnomad')
#' ACMG.exac <- NULL; ACMG.gnomad <- NULL
#' for (gene in ACMG.panel) {
#'   print(sprintf("[%d/%d] %s",which(gene==ACMG.panel),length(ACMG.panel),gene))
#'   ACMG.exac <- rbind(ACMG.exac,import_file_exac(gene, "exac"))
#'   ACMG.gnomad <- rbind(ACMG.gnomad,import_file_exac(gene, "gnomad"))
#' }
#' ACMG.exac <- ACMG.exac[!duplicated(ACMG.exac$VAR_ID),]
#' ACMG.gnomad <- ACMG.gnomad[!duplicated(ACMG.gnomad$VAR_ID),]
#' }

import_file_exac <- function(gene, dataset) {
  file_name <- sprintf("%s/%s_%s.csv", dataset, dataset, gene)
  output <- read.csv(file_name, stringsAsFactors = FALSE)
  output$Number.of.Hemizygotes <- NULL #Inconsistently present column; removal allows row aggregation
  # Correcting for some alternate naming conventions
  if ("Conseq." %in% colnames(output))
    output <- output %>% rename(Consequence = Conseq.)
  if ("Count" %in% colnames(output))
    output <- output %>% rename(Allele.Count = Count)
  # Imputing missing South Asian values for NF2
  if (!("Allele.Count.South.Asian" %in% colnames(output))) {
    output$Allele.Number.South.Asian <- (2*output$Allele.Number) -
      (output %>% select(contains("Allele.Number")) %>% rowSums)
    output$Allele.Count.South.Asian <- (2*output$Allele.Count) - 
      (output %>% select(contains("Allele.Count")) %>% rowSums)
    output$Homozygote.Count.South.Asian <- (2*output$Number.of.Homozygotes) - 
      (output %>% select(contains("Homozygote")) %>% rowSums)
  }
  output <- cbind(GENE = gene, output[nchar(paste(output$Alternate,output$Reference))==3,]) %>% 
    select(GENE, AF_EXAC = contains("Freq"), CHROM=Chrom, POS=Position, 
           ID=RSID, REF=Reference, ALT=Alternate, Annotation = contains("Annot"), everything()) %>% 
    unite(VAR_ID, CHROM, POS, REF, ALT, sep = "_", remove = F) %>% arrange(VAR_ID)
  tags <- list("African","Latino","East.Asian","European","South.Asian")
  european <- output %>% select(contains("Finnish"), contains("European"))
  if (dataset == "gnomad") {
    european <- output %>% select(contains("Finnish"), contains("European"), contains("Jewish"))
    output <- output %>% select(GENE, AF_GNOMAD = AF_EXAC, everything())
  }
  output$Allele.Count.European <- european %>% select(contains("Allele.Count")) %>% rowSums
  output$Allele.Number.European <- european %>% select(contains("Allele.Number")) %>% rowSums
  exac_af <- output[,sprintf("Allele.Count.%s", tags)] / output[,sprintf("Allele.Number.%s", tags)]
  colnames(exac_af) <- sprintf("AF_%s_%s", toupper(dataset), c("AFR","AMR","EAS","EUR","SAS"))
  output <- cbind(output, exac_af) %>% 
    select(GENE, contains(toupper(dataset)), everything())
  if (rm_rs1805124)
    output <- output[output$ID!="rs1805124",]
  return(output)
}