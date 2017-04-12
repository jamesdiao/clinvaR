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
#' @export

import_file_exac <- function(gene, dataset) {
  dir <- system.file("extdata", "", package = "clinvaR")
  file_name <- sprintf("%s%s/%s_%s.csv", dir, dataset, dataset, gene)
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
  #if (rm_rs1805124)
    output <- output[output$ID!="rs1805124",]
  return(output)
}



#' Imports downloaded VCF from 1000 Genomes
#'
#' This function imports a gene's VCF and extracts attributes from INFO 
#' into a specified format for other clinvaR analyses
#' @usage import_file_1000g(gene)
#' @examples import_file_1000g('BRCA2')
#' ACMG.1000g <- NULL
#' for (gene in ACMG.panel) {
#' ACMG.1000g <- rbind(ACMG.1000g,import_file_1000g(gene))
#' }
#' #ACMG.1000g[duplicated(ACMG.1000g$VAR_ID),1:8] # Display and remove duplicates
#' ACMG.1000g <- ACMG.1000g[!duplicated(ACMG.1000g$VAR_ID),]
#' @export

import_file_1000g <- function(gene) {
  map <- system.file("extdata", "Supplementary_Files/phase3map.txt", package = "clinvaR") %>% 
    read.table(stringsAsFactors = F, header = T) %>% as.data.frame
  header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", 
              "FILTER", "INFO", "FORMAT", as.character(map$sample))
  dir <- system.file("extdata", "", package = "clinvaR")
  output <- read.table(sprintf('%s1000G/%s_genotypes.vcf', dir, gene), stringsAsFactors = FALSE)
  #Add header
  names(output)[1:length(header)] <- header
  #Remove all single alt indels
  output <- output[nchar(output$REF)==1,] #deletions
  alt_num <- sapply(strsplit(output$ALT,","),length) #number of alts
  acceptable_nchar <- 2*alt_num-1 #adds in the length from commas, if each alt is 1 nt.
  output <- output[nchar(output$ALT)==acceptable_nchar,] #insertions
  alt_num <- sapply(strsplit(output$ALT,","),length) #recalculate
  paired = which(alt_num!=1) #all with ,
  #Add AF Column
  af <- strsplit(output$INFO,";") %>% sapply("[", 2) %>% 
    strsplit("AF=") %>% sapply("[", 2) %>% strsplit(",") %>% sapply(as.numeric)
  output <- cbind(GENE = gene, "AF_1000G"=I(af), output) #Places it at the front of output
  front_cols <- 1:(grep("HG00096",colnames(output))-1)
  if (length(paired)!=0) {
    #Limit max vector length by sapply(strsplit(output$ALT,","),length)
    sapply(paired, function(rownum) { #For every row
      sapply(as.character(1:alt_num[rownum]), function(num) {
        grepl(paste(num,"|",sep = ""), output[rownum,-front_cols], fixed=T) +
          grepl(paste("|",num,sep = ""), output[rownum,-front_cols], fixed=T)
      }) %>% t -> temp
      split(temp, rep(1:ncol(temp), each = nrow(temp))) %>% setNames(NULL) 
      #Separate into list of vectors (1 entry for counting each ALT)
    }) %>% t -> insert
    insert <- cbind(output[paired,front_cols],insert)
    colnames(insert) <- colnames(output)
    insert <- insert %>% #adds front_col info
      mutate(ALT = strsplit(ALT,",")) %>% #Splits ALTS
      unnest() %>% #Unnests everything
      select(GENE, AF_1000G, CHROM, POS, ID, REF, ALT, everything()) #Reorders everything
    output <- output[-paired,] #Removes paired
  }
  output <- cbind(output[,front_cols],
                  apply(output[,-front_cols], 2, function(y) {
                    grepl("1|", y, fixed=T) +
                      grepl("|1", y, fixed=T)
                  }) ) #convert to logical
  if (length(paired)!=0)
    output <- rbind(output, insert) #joins the two
  output$AF_1000G <- as.numeric(output$AF_1000G)
  unite(output, VAR_ID, CHROM, POS, REF, ALT, sep = "_", remove = F) %>% arrange(VAR_ID)
  #Make VAR_ID, arrange by VAR_ID
}