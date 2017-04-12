#' Downloads VCF from 1000 Genomes using tabix
#'
#' This function downloads a gene's VCF from the 1000 Genomes FTP based on its transcription region
#' as listed in the UCSC Genome database. 
#' @usage download_1000g(gene, download)
#' @examples download_1000g('BRCA2', download = T)
#' download_output <- sapply(ACMG.panel, function(gene) download_1000g(gene, download = T)) %>% t
#' @export


download_1000g <- function(gene, download) {
  if (missing(download)) {
    download <- TRUE
  }
  #for tracking: #gene %>% paste(which(ACMG.panel==gene)) %>% paste(length(ACMG.panel), sep = "/") %>% print
  success <- FALSE
  refGene <- sprintf("select * from refGene where name2 = \"%s\" limit 20", gene) %>% query
  UCSC <- select(refGene, name, chrom, start = txStart, end = txEnd)
  if (nrow(UCSC) == 0) { #No hit on refGene
    return(rep("NOT_FOUND",5) %>% setNames(c("name","chrom","start","end","downloaded")))
  } else {
    if (nrow(UCSC) > 1) #Multiple hits: take the widest range
      UCSC <- UCSC[which.max(UCSC$end-UCSC$start),]
    if (download) {
      # gets [n] from chr[n]
      chrom.num <- strsplit(UCSC$chrom, split = "chr")[[1]][2]
      # different version for chromosomes X and Y
      version <- switch(chrom.num, "X" = "shapeit2_mvncall_integrated_v1b",
                        "Y" = "integrated_v2a", "shapeit2_mvncall_integrated_v5a")
      command <- paste("tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.%s.",
                       "phase3_%s.20130502.genotypes.vcf.gz %s:%s-%s > %s_genotypes.vcf", sep = "")
      sprintf(command, UCSC$chrom, version, chrom.num, UCSC$start, UCSC$end, gene) %>% system
      Sys.sleep(2)
      # Checks whether the file exists and has non-zero size
      exists <- grepl(paste(gene,"_genotypes.vcf",sep =""), system("ls", intern = T)) %>% sum > 0
      file.size <- strsplit(paste("stat ","_genotypes.vcf", sep = gene) %>% 
                              system(intern = T), " ")[[1]][8]
      success <- exists & file.size > 0
    }
  }
  return(c(UCSC,"downloaded" = success))
}