#' Downloads VCF from 1000 Genomes using tabix
#'
#' This function downloads a gene's VCF from the 1000 Genomes FTP based on its transcription region
#' as listed in the UCSC Genome database.
#' @usage download_1000g(genes, download)
#' @param genes character; name of a gene(s) indexed in refGene
#' @param download logical; indicates whether to download the file. 
#' FALSE still returns transcription coordinates and RefGene details from UCSC. 
#' @examples 
#' download_1000g('BRCA2', download = T)
#' download_output <- sapply(ACMG.panel, function(gene) download_1000g(gene, download = T)) %>% t
#' @export

download_1000g <- function(genes, download, verbose) {
  if (missing(verbose)) 
    verbose <- TRUE
  if (missing(download)) 
    download <- TRUE 
  dir <- system.file("extdata", "1000G", package = "clinvaR")
  refGene <- readRDS(sprintf('%s/refGene.rds', dir))
  ptm <- proc.time()
  download_output <- sapply(genes, function(gene) {
    success <- FALSE
    elapsed <- h_read((proc.time() - ptm)['elapsed'])
    if (length(genes) > 1 & verbose) {
      print(sprintf("Downloading [%d/%d] %s, Time Elapsed: %s", 
                    which(gene==genes), length(genes), gene, elapsed), 
            quote = F)
    } else {
      print(sprintf("Downloading %s, Time Elapsed: %s", gene, elapsed), 
            quote = F)
    }
    
    if (gene %in% refGene$gene) { 
      UCSC <- refGene[refGene$gene == gene, -1]
    } else {
      return(rep("NOT_FOUND",5) %>% setNames(c("name","chrom","start","end","downloaded")))
    }
    
    if (download) {
      # different version for chromosomes X and Y
      chrom <- UCSC$chrom %>% replace(UCSC$chrom == "23", "X") %>% replace(UCSC$chrom == "24", "Y")
      version <- switch(as.character(chrom), "X" = "shapeit2_mvncall_integrated_v1b",
                        "Y" = "integrated_v2a", "shapeit2_mvncall_integrated_v5a")
      file.name <- sprintf(
                      fmt = paste0("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/",
                             "ALL.%s.phase3_%s.20130502.genotypes.vcf.gz"), 
                      paste0('chr', chrom), version
                   )
      output <- tabix.read(tabixFile = file.name, 
                           tabixRange = sprintf('%s:%s-%s', chrom, UCSC$start, UCSC$end)
                           )
      system("rm *.genotypes.vcf.gz.tbi")
      output <- read.table(text = paste(output, collapse = "\n"), header = F, stringsAsFactors = F, 
                          comment.char = "", quote = "", sep = "\t")
      saveRDS(output, sprintf('%s/%s_genotypes_vcf.rds', dir, gene))
      
      # Checks whether the file exists and has non-zero size
      exists <- grepl(sprintf("%s_genotypes_vcf.rds", gene), 
                      system(sprintf("ls %s", dir), intern = T)
                      ) %>% any #try file.exists
      path <- sprintf("%s/%s_genotypes_vcf.rds", dir, gene)
      success <- file.exists(path) & (file.size(path) > 0)
    }
    return(c(UCSC,"downloaded" = success))
  })
  return(download_output)
}

#' Delete All Downloaded VCF.rds Files
#'
#' This function deletes all files generated using download_1000g()
#' 
#' @usage delete_downloads()
#' @export
#' 

delete_downloads <- function() {
  dir <- system.file("extdata", package = "clinvaR")
  contents <- system(sprintf('ls %s/1000G', dir), intern = T)
  contents <- contents[grepl('_genotypes_vcf.rds', contents)]
  file.remove(contents)
}