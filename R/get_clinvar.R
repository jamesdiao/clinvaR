#' Import ClinVar VCF
#'
#' This function allows you to read a ClinVar VCF 
#' into a table and extract important information
#' from the INFO section. 
#' @usage get_clinvar(clinvar_file)
#' @param clinvar_file Specifies the path to the VCF. 
#' @examples 
#' clinvar_date <- system("ls ClinVar_Reports/VCF/clinvar*.vcf.gz", intern = T) %>% 
#' strsplit("clinvar_") %>% sapply("[[",2) %>% unlist %>% 
#' strsplit(".vcf.gz")  %>% sapply("[[",1) %>% unlist %>% as.Date() %>% max
#' # collects the most recent date in './ClinVar_Reports/VCF'
#' clinvar_file <- sprintf("ClinVar_Reports/VCF/clinvar_%s.vcf.gz", clinvar_date)
#' # find the file name
#' clinvar <- get_clinvar(clinvar_file)
#' clinvar <- clinvar[!duplicated(clinvar$VAR_ID),] # remove duplicates
#' @export

get_clinvar <- function(clinvar_file) {
  file.by.line <- readLines(clinvar_file)
  #file_date <- as.Date(strsplit(file.by.line[2],"=")[[1]][2], "%Y%m%d")
  #system(sprintf("mv %s ClinVar_Reports/clinvar_%s.vcf", clinvar_file, file_date))
  clean.lines <- file.by.line[!grepl("##.*", file.by.line)] #Remove ## comments
  clean.lines[1] <- sub('.', '', clean.lines[1]) #Remove # from header
  input <- read.table(text = paste(clean.lines, collapse = "\n"), header = T, stringsAsFactors = F, 
                      comment.char = "", quote = "", sep = "\t")
  input <- input[nchar(input$REF)==1,] #deletions
  alt_num <- sapply(strsplit(input$ALT,","),length) #number of alts
  acceptable_nchar <- 2*alt_num-1 #adds in the length from commas, if each alt is 1 nt.
  input <- input[nchar(input$ALT)==acceptable_nchar,] #insertions
  input$ALT <- strsplit(input$ALT,",")
  split_all <- strsplit(input$INFO,";")
  has.clndsdb <- any(grepl("CLNDSDB", split_all))
  has.clndsdbid <- any(grepl("CLNDSDBID", split_all))
  has.clnrevstat <- any(grepl("CLNREVSTAT", split_all))
  
  split_info <- function(name) {
    sapply(split_all, function(entry) {
      entry[grep(name,entry)]
    }) %>% strsplit("=") %>% sapply(function(x) x[2]) %>% strsplit(",")
  }
  input$CLNALLE <- split_info("CLNALLE=") %>% sapply(as.integer)
  input$CLNSIG <- split_info("CLNSIG=")
  input$CLNDBN <- split_info("CLNDBN=")
  if (has.clnrevstat)
    input$CLNREVSTAT <- split_info("CLNREVSTAT=")
  if (has.clndsdb)
    input$CLNDSDB <- split_info("CLNDSDB=")
  if (has.clndsdbid)
    input$CLNDSDBID <- split_info("CLNDSDBID=")
  #CLNALLE has 0,-1,3,4 --> CLNSIG has 1,2,3,4 --> ALT has 1. 
  taking <- sapply(input$CLNALLE, function(x) x[x>0] ) #Actual elements > 0. Keep these in CLNSIG and ALT 
  taking_loc <- sapply(input$CLNALLE, function(x) which(x>0) )#Tracks locations for keeping in CLNALLE
  keep <- sapply(taking, length)>0 #reduce everything to get rid of 0 and -1
  # Reduce, reduce, reduce. 
  taking <- taking[keep]
  taking_loc <- taking_loc[keep]
  input <- input[keep,]
  
  #Make this more readable
  input$ALT <- sapply(1:nrow(input), function(row) {
    input$ALT[[row]][taking[[row]]]
  })
  
  col_subset <- function(name) {
    sapply(1:nrow(input), function(row) {
      unlist(input[row,name])[taking_loc[[row]]]
    })
  }
  input$CLNSIG <- col_subset("CLNSIG")
  input$CLNALLE <- col_subset("CLNALLE")
  input$CLNDBN <- col_subset("CLNDBN")
  if (has.clnrevstat)
    input$CLNREVSTAT <- col_subset("CLNREVSTAT")
  if (has.clndsdb)
    input$CLNDSDB <- col_subset("CLNDSDB")
  if (has.clndsdbid)
    input$CLNDSDBID <- col_subset("CLNDSDBID")
  filter_condition <- input[,unlist(lapply(input, typeof))=="list"] %>% 
    apply(1,function(row) !any(is.na(row)))
  input <- input %>% filter(filter_condition) %>%
    unnest %>% unite(VAR_ID, CHROM, POS, REF, ALT, sep = "_", remove = F) %>%
    select(VAR_ID, CHROM, POS, ID, REF, ALT, CLNALLE, CLNSIG, everything()) %>% 
    mutate(CLNSIG = strsplit(CLNSIG,"|",fixed = T)) %>% 
    mutate(CLNDBN = strsplit(CLNDBN,"|",fixed = T)) %>% 
    mutate(POS = as.integer(POS))
  if (has.clnrevstat)
    input <- input %>% mutate(CLNREVSTAT = strsplit(CLNREVSTAT,"|",fixed = T))
  if (has.clndsdb)
    input <- input %>% mutate(CLNDSDB = strsplit(CLNDSDB,"|",fixed = T)) 
  if (has.clndsdbid)
    input <- input %>% mutate(CLNDSDBID = strsplit(CLNDSDBID,"|",fixed = T)) 
  input$CLNSIG <- sapply(input$CLNSIG, function(x) as.integer(x))
  input$INTERP <- sapply(input$CLNSIG, function(x) any(x %in% c(4,5)) & !(any(x %in% c(2,3)))) 
  input
}