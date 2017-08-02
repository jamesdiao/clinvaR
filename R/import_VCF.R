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
  
  if (missing("genes")) {
    contents <- system(sprintf('ls %s/1000G', dir), intern = T)
    contents <- contents[grepl('_genotypes_vcf.rds', contents)]
    if (length(contents)==0) {
      print('Error: no genes detected in downloads folder', quote = F)
      return(NULL)
    } else {
      genes <- str_match(string = contents, pattern = '([^_]*)_genotypes_vcf.rds')[,2]
    }
  } else {
    genes <- toupper(genes)
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
      print(sprintf("Warning at %s: VCF file not found", 
                    which(gene==genes), length(genes), gene), quote = F)
      exists <- download_1000g(gene)['downloaded',] %>% unlist
      print(ifelse(exists, "SUCCESS", "FAILURE"))
    } 
    if (exists) {
      print(sprintf("Importing [%d/%d] %s", 
                    which(gene==genes), length(genes), gene), quote = F)
      combined <- rbind(combined, temp_function(gene))
    }
  }
  combined <- combined[!duplicated(combined$VAR_ID),]
  class(combined) <- c("plain_vcf", class(combined))
  return(combined)
}

#' Compute and Plot Aggregate Allele Frequencies for 1000 Genomes
#'
#' This function takes a plain ClinVar-sequencing dataset and 
#' returns an ancestry-stratified plot of counts or frequencies.    
#' @usage var_plot_1000g(vcf)
#' var_plot_1000g(vcf, fraction = FALSE)
#' @param vcf plain_vcf data.frame; clinvaR-processed VCF from 1000 Genomes. 
#' @param fraction logical; if TRUE, plots fraction with a finding. If FALSE, plots total counts. 
#' @export

plot.plain_vcf <- function(vcf, fraction, sd) {
  if (missing(fraction)) 
    fraction <- FALSE
  if (missing(sd))
    sd <- TRUE
  if (!is.logical(fraction)) {
    print("Fraction must be logical: set as FALSE")
    fraction <- FALSE
  }
  return(var_plot_1000g(vcf, fraction, sd))
}
