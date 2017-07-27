## clinvaR
R Package for for gene-level analysis of variants in the ClinVar database. 

#### Purpose
The goal of `clinvaR` is to quickly generate a VCF of annotated variants relevant to a list of genes, and produce simple, ancestry-informed analyses. 

#### Installation
``` r
# install.packages("devtools")
devtools::install_github("jamesdiao/clinvaR")
require(clinvaR)
```

#### Usage
To acquire a VCF with the July 5, 2017 ClinVar file and all default settings, simply run `annotate(genes)`.  

More specific modifications can be made by changing the 4 steps:

1.  Select a list of genes that you are interested in. This can be done by manual input, by importing a .tsv from MacArthur's gene lists (stored locally), or by importing your own .tsv file: `get_genes()`
2.  Download and import relevant variants from 1000 Genomes: `download_1000g(), import_file_1000g()`. 
3.  Select a version of ClinVar. These can be downloaded directly: `download_clinvar()`. Alternatively, ClinVar VCFs before July 2017 can be retrieved locally: `get_clinvar()`.
4.  Use the ClinVar file to annotate the gene-level variant VCF: `annotate_1000g()`.


`clinvaR` also provides a few basic analysis and visualization functions that can be run directly : `freq_over_time`, `var_plot_1000g`. 

###### James Diao | Arjun Manrai | Kohane Lab | 07/27/2017