% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clinvaR.R
\docType{package}
\name{clinvaR-package}
\alias{clinvaR}
\alias{clinvaR-package}
\title{clinvaR: collection, analysis, and visualization tools for ClinVar data}
\description{
clinvaR rapidly generates a VCF of annotated variants for a list of genes and 
produces simple, ancestry-informed analyses.
}
\details{
To acquire a VCF with the July 5, 2017 ClinVar file and all default settings, 
simply run annotate(genes). More specific modifications can be made by changing the 4 steps:
\itemize{
\item  Select a list of genes that you are interested in. This can be done by manual input, 
by importing a .tsv from MacArthur's gene lists (stored locally), or by importing your own 
.tsv file: get_genes()
\item Download and import relevant variants from 1000 Genomes: 
download_1000g(), import_file_1000g().
\item Select a version of ClinVar. These can be downloaded directly: download_clinvar(). 
Alternatively, ClinVar VCFs before July 2017 can be retrieved locally: get_clinvar().
\item Use the ClinVar file to annotate the gene-level variant VCF: annotate_1000g().
}

clinvaR also provides a few basic analysis and visualization functions that can be run directly : freq_over_time, var_plot_1000g.
To learn more about clinvaR, start with the vignettes:
`browseVignettes(package = "clinvaR")`
}
\author{
\strong{Maintainer}: James Diao \email{diao.james@gmail.com}

Other contributors:
\itemize{
  \item Arjun Manrai \email{manrai@post.harvard.edu} [thesis advisor]
}

}
