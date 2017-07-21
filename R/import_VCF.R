#' Imports downloaded VCF from 1000 Genomes
#'
#' This function imports a gene's VCF and extracts attributes from INFO 
#' into a specified format for other clinvaR analyses
#' @usage import_file_1000g(genes, dataset, path)
#' @param gene character; any of the pre-stored genes in extdata/1000G
#' Not case-sensitive. Defaults to 'gnomAD'.
#' @param path character; path to an importable VCF. Overrides gene input. 
#' @examples import_file_1000g('BRCA2')
#' ACMG.1000g <- NULL
#' for (gene in ACMG.panel) {
#' ACMG.1000g <- rbind(ACMG.1000g,import_file_1000g(gene))
#' }
#' #ACMG.1000g[duplicated(ACMG.1000g$VAR_ID),1:8] # Display and remove duplicates
#' ACMG.1000g <- ACMG.1000g[!duplicated(ACMG.1000g$VAR_ID),]
#' @export

import_file_1000g <- function(genes) {
  dir <- system.file("extdata", package = "clinvaR")
  
  if (missing(genes)) {
    contents <- system(sprintf('ls %s/1000G', dir), intern = T)
    contents <- contents[grepl('_genotypes_vcf.rds', contents)]
    if (length(contents)==0) {
      stop('Error: no genes detected in downloads folder')
    } else {
      genes <- str_match(string = contents, pattern = '([^_]*)_genotypes_vcf.rds')[,2]
    }
  }
  
  temp_function <- function(gene) {
    path <- sprintf("%s/1000G/%s_genotypes_vcf.rds", dir, gene)
    output <- readRDS(path)
    colnames(output) <- header
    
    #Remove all single alt indels
    output <- output %>% filter(nchar(REF)==1) %>% #deletions
      filter(nchar(ALT)==2*str_count(ALT,",")+1) #insertions
    alt_num <- str_count(output$ALT,",") + 1 #number of alternates
    paired = which(alt_num!=1) #all with ,
    #Add AF Column
    af <- str_match(output$INFO, ';AF=([^;]*);')[,2] %>% strsplit(",") %>% sapply(as.numeric)
    output <- data.frame(GENE = gene, "AF_1000G"=I(af), output) #Places it at the front of output
    front_cols <- 1:(grep("HG00096",colnames(output))-1)
    if (length(paired)!=0) {
      #Limit max vector length by sapply(strsplit(output$ALT,","),length)
      sapply(paired, function(rownum) { #For every row
        sapply(as.character(1:alt_num[rownum]), function(num) {
          grepl(paste0(num, "|"), output[rownum, -front_cols], fixed=T) +
            grepl(paste0("|", num), output[rownum, -front_cols], fixed=T)
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
    output <-  data.frame(output[,front_cols],
                          apply(output[,-front_cols], 2, function(y) {
                            grepl("1|", y, fixed=T) +
                              grepl("|1", y, fixed=T)
                          }) ) #convert to logical
    if (length(paired)!=0)
      output <- rbind(output, insert) #joins the two
    output <- unite(output, VAR_ID, CHROM, POS, REF, ALT, sep = "_", remove = F) %>% 
      arrange(VAR_ID) %>% mutate(AF_1000G = AF_1000G %>% unlist)
    return(output)
  }
  
  map <- read.table(file = sprintf("%s/Supplementary_Files/phase3map.txt", dir), stringsAsFactors = F, header = T) %>% as.data.frame
  header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", as.character(map$sample))
  combined <- NULL
  for (gene in genes) {
    exists <- any(grepl(gene, system(sprintf('ls %s/1000G', dir), intern = T)))
    if (!exists) {
      print(sprintf("ERROR at [%d/%d] %s: VCF file not found", which(gene==genes), length(genes), gene))
      exists <- download_1000g(gene)['downloaded',] %>% unlist
      print(ifelse(exists, "SUCCESS", "FAILURE"))
    } 
    if (exists) {
      print(sprintf("Importing [%d/%d] %s", which(gene==genes), length(genes), gene))
      combined <- rbind(combined, temp_function(gene))
    }
  }
  combined <- combined[!duplicated(combined$VAR_ID),]
  return(combined)
}

#' Imports downloaded VCF from ExAC or gnomAD
#'
#' This function imports a gene's VCF and extracts attributes into a 
#' specified format for other clinvaR analyses
#' 
#' @usage import_file_exac(gene, dataset, path)
#' @param gene character; any of the pre-stored genes in extdata/ExAC or extdata/gnomAD
#' @param dataset character; choose from 'ExAC' or 'gnomAD'. 
#' Not case-sensitive. Defaults to 'gnomAD'.
#' @param path character; path to an importable VCF. Overrides gene and dataset inputs. 
#' @examples 
#' import_file_exac('BRCA2', 'gnomad')
#' ACMG.exac <- NULL; ACMG.gnomad <- NULL
#' for (gene in ACMG.panel) {
#'   print(sprintf("[%d/%d] %s",which(gene==ACMG.panel),length(ACMG.panel),gene))
#'   ACMG.exac <- rbind(ACMG.exac,import_file_exac(gene, "exac"))
#'   ACMG.gnomad <- rbind(ACMG.gnomad,import_file_exac(gene, "gnomad"))
#' }
#' ACMG.exac <- ACMG.exac[!duplicated(ACMG.exac$VAR_ID),]
#' ACMG.gnomad <- ACMG.gnomad[!duplicated(ACMG.gnomad$VAR_ID),]
#' #@export

import_file_exac <- function(gene, dataset, path) {
  dir <- system.file("extdata", "", package = "clinvaR")
  if (missing(dataset)) dataset <- 'gnomAD'
  if (missing(path)) {
    output <- sprintf("%s%s/%s_%s.csv", dir, dataset, dataset, gene) %>% 
      read.csv(stringsAsFactors = FALSE)
  } else {
    output <- read.csv(path, stringsAsFactors = FALSE)
  }
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
