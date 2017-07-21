#' Imports List of Genes from File
#'
#' This function imports a list of genes from a specified file location, or from the MacArthur gene lists. 
#' @usage get_genes(file)
#' @param file character; specified file location, or local file names from extdata/MacArthur_Gene_Lists/gene_lists/lists
#' @examples get_genes('lmm_hcm.tsv')
#' @export

get_genes <- function(file) {
  if (file.exists(file)) {
    return(read.table(file) %>% unlist() %>% unique() %>% as.character())
  } else {
    if (!grepl('.tsv', file, fixed = T))
      file <- paste0(file, '.tsv')
    file <- system.file(sprintf('extdata/MacArthur_Gene_Lists/lists/%s', file), 
                        package = 'clinvaR') %>% tolower()
    if (nchar(file) > 0) {
      return(read.table(file) %>% unlist() %>% unique() %>% as.character())
    } else {
      stop('ERROR: file not found.')
    }
    
  }
}
