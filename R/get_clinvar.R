#' Import ClinVar VCF
#'
#' This function allows you to read a ClinVar VCF into a table 
#' and extract important information from the INFO section. 
#' 
#' @usage download_clinvar()
#' download_clinvar(file)
#' @param file character; specifies the path to the VCF. 
#' @details getclinvar() collects the latest VCF from 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz'
#' Specifying a file path will collect the VCF from that path instead. 
#' @examples 
#' clinvar <- sprintf("ClinVar_Reports/VCF/clinvar_%s.vcf.gz", clinvar_date) %>% 
#'               download_clinvar(file)
#' clinvar <- clinvar[!duplicated(clinvar$VAR_ID),] # remove duplicates
#' @export

download_clinvar <- function(file) {
  
  newest_clinvar <- function() {
    dir <- system.file("extdata", "Supplementary_Files/", package = "clinvaR")
    file <- sprintf("%sclinvar_%s.vcf.gz", dir, Sys.Date())
    download.file(url = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz", 
                  destfile = file, method = "internal")
    system(sprintf("gunzip %s", file))
    file <- sprintf("%sclinvar_%s.vcf", dir, Sys.Date())
    return(file)
  }
  
  if (missing(file))
    file <- newest_clinvar()
  if (!file.exists(file))
    file <- newest_clinvar()
  
  extract_element <- function(phrase) {
    str_match_all(input$INFO, sprintf('%s=([^;]*);', phrase)) %>% 
      lapply('[[', 2) %>% unlist
  }
  file.by.line <- readLines(file)
  #file_date <- as.Date(strsplit(file.by.line[2],"=")[[1]][2], "%Y%m%d")
  #system(sprintf("mv %s ClinVar_Reports/clinvar_%s.vcf", file, file_date))
  clean.lines <- file.by.line[!grepl("##.*", file.by.line)] #Remove ## comments
  clean.lines[1] <- sub('.', '', clean.lines[1]) #Remove # from header
  input <- read.table(text = paste(clean.lines, collapse = "\n"), header = T, stringsAsFactors = F, 
                      comment.char = "", quote = "", sep = "\t")
  input <- input[nchar(input$REF)==1,] #deletions
  alt_num <- sapply(strsplit(input$ALT,","),length) #number of alts
  acceptable_nchar <- 2*alt_num-1 #adds in the length from commas, if each alt is 1 nt.
  input <- input[nchar(input$ALT)==acceptable_nchar,] #insertions
  input$ALT <- strsplit(input$ALT,",")
  input$CLNALLE <- extract_element('CLNALLE') %>% strsplit(',', fixed = T) %>% lapply(as.integer)
  input$CLNSIG <- extract_element('CLNSIG') %>% strsplit(',', fixed = T)
  input$CLNDBN <- extract_element('CLNDBN') %>% strsplit(',', fixed = T)
  
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
  filter_condition <- input[,unlist(lapply(input, typeof))=="list"] %>% 
    apply(1,function(row) !any(is.na(row)))
  input <- input %>% filter(filter_condition) %>%
    unnest %>% unite(VAR_ID, CHROM, POS, REF, ALT, sep = "_", remove = F) %>%
    select(VAR_ID, CHROM, POS, REF, ALT, ID, CLNSIG, CLNDBN) %>% 
    mutate(CLNSIG = strsplit(CLNSIG,"|",fixed = T)) %>% 
    mutate(CLNDBN = strsplit(CLNDBN,"|",fixed = T)) %>% 
    mutate(POS = as.integer(POS))
  input$CLNSIG <- sapply(input$CLNSIG, function(x) as.integer(x))
  input$pathogenic_incl_conflicts <- sapply(input$CLNSIG, function(x) any(x %in% c(4,5)))
  input$pathogenic_no_conflicts <- sapply(input$CLNSIG, function(x) any(x %in% c(4,5)) & !(any(x %in% c(2,3)))) 
  #input$INTERP <- input$pathogenic_no_conflicts
  #input$LMM <- grepl("Laboratory_for_Molecular_Medicine",input$INFO)
  return(input)
}

#' Get List of Available ClinVar Version Dates
#'
#' This function returns the dates of all available ClinVar VCFs
#' 
#' @usage get_date_list() 
#' @export

get_date_list <- function(resolution) {
  if (missing(resolution)) 
    resolution <- "all"
  dir <- system.file("extdata", package = "clinvaR")
  clinvar_reports <- system(sprintf("ls %s/Archive_Tables", dir), intern = T)
  clinvar_reports <- clinvar_reports[grep(".tsv",clinvar_reports)]
  do.call("rbind", lapply(clinvar_reports, function(tsv) {
    read.table(file = sprintf("%s/Archive_Tables/%s", dir, tsv), sep = "\t", 
               header = T, stringsAsFactors = F)
  })) -> archive
  all_dates <- regmatches(archive$Name,regexpr("^File:clinvar_(20.{6})\\.vcf\\.gz$",archive$Name)) %>% 
    str_extract("20.{6}") %>% as.Date(format = "%Y%m%d") %>% unique()
  if (resolution == "year")
    return(get_closest_date(sprintf("%s-01-01", 2012:2017)))
  return(all_dates)
}

#' Get Closest ClinVar Version Dates
#'
#' This function returns the closest ClinVar date to the input date. 
#' If there is no input, it defaults to the current date. 
#' 
#' @usage get_date_list() 
#' get_date_list(dates)
#' @param file character; input date in the format \%y-\%m-\%d. Ex: 2017-07-05.  
#' @export

get_closest_date <- function(dates) {
  if (missing(dates)) 
    dates <- Sys.Date()
  all_dates <- get_date_list()
  sapply(dates, function(date) {
    as.character(all_dates)[which.min(abs(all_dates - as.Date(date)))]
  })
}

#' Get ClinVar Version from Date
#'
#' This function returns the ClinVar version from the date closest to the input date. 
#' If there is no input, it defaults to the current date. 
#' 
#' @usage get_clinvar()
#' get_clinvar(date)
#' @param file character; input date in the format \%y-\%m-\%d. Ex: 2017-07-05.  
#' @export

get_clinvar <- function(date) {
  if (missing(date)) 
    date <- Sys.Date()
  system.file("extdata", sprintf('RDS/clinvar_%s.rds', get_closest_date(date)), 
              package = "clinvaR") %>% readRDS %>% return
}

  
  