#' Imports downloaded VCF from 1000 Genomes
#'
#' This function imports a gene's VCF and extracts attributes from INFO 
#' into a specified format for other clinvaR analyses
#' @usage import_file_1000g(gene)
#' @examples import_file_1000g('BRCA2')
#' ACMG.1000g <- NULL
#' header <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", as.character(map$sample))
#' for (gene in ACMG.panel) {
#' ACMG.1000g <- rbind(ACMG.1000g,import_file_1000g(gene))
#' }
#' #ACMG.1000g[duplicated(ACMG.1000g$VAR_ID),1:8] # Display and remove duplicates
#' ACMG.1000g <- ACMG.1000g[!duplicated(ACMG.1000g$VAR_ID),]

import_file_1000g <- function(gene) {
  #for tracking: 
  sprintf("%s [%s/%s]", gene, grep(gene, ACMG.panel), length(ACMG.panel)) %>% print(quote = F)
  name <- paste("1000G",paste(gene,"genotypes.vcf", sep = "_"), sep = "/")
  output <- read.table(paste(getwd(),name,sep="/"), stringsAsFactors = FALSE)
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