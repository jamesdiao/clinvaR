#' Connect to UCSC hg19 Genome Database
#'
#' This function makes a connection to the hg19 genome database and defines a query function. 
#' @usage connect_UCSC
#' @export

connect_UCSC <- function() {
  for (con in dbListConnections(MySQL())) dbDisconnect(con)
  con <- dbConnect(MySQL(), user = 'genome',
                   dbname = 'hg19', host = 'genome-mysql.cse.ucsc.edu',
                   unix.sock = "/Applications/MAMP/tmp/mysql/mysql.sock")
  query <- function (input) { suppressWarnings(dbGetQuery(con, input)) }
}