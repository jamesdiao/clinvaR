This is a vignette intended to introduce new users to `clinvaR`
(<https://github.com/jamesdiao/clinvaR>).

Purpose
-------

For `clinvaR`, the first goal is to go from genes to a vcf of related
variants, along with classifications from any version of ClinVar.  
This involves 3 steps:

1.  Select a list of genes that you are interested in. clinvaR makes it
    easy to download and import relevant variants from 1KG.  
2.  Select a version of ClinVar to assign pathogenicity labels with. All
    previous ClinVar VCFs are stored as binary files and can be
    retrieved with a "closest date" function.  
3.  Select an analysis method. This part is not ready yet.

Dependencies
------------

clinvaR has a number of file dependencies that are generally useful for
these types of projects. I have attached them here.

    ## [1] remotes pander  ggplot2 tibble  tidyr   dplyr   stringr

Select a List of Genes
----------------------

The first step is to have a list of genes that you are interested in.

### Method 1: Direct assignment

The easiest way is to directly save your genes of interest in
gene\_panel as a character vector (all caps). We have taken the example
of the 20 genes on the Laboratory of Molecular Medicine's gene testing
panel for hypertrophic cardiomyopathy.

    gene_panel <- c("ACTC1", "ACTN2", "CSRP3", "GLA", "LAMP2", "MYBPC3", "MYH7", 
    "MYL2", "MYL3", "MYOZ2", "NEXN", "PLN", "PRKAG2", "PTPN11", "RAF1", 
    "TNNC1", "TNNI3", "TNNT2", "TPM1", "TTR")

### Method 2: Import from text file

#### From own file

Another way is to save your gene list in a text file (1 gene per line)
and import it.

    path <- system.file("extdata/MacArthur_Gene_Lists/lists/lmm_hcm.tsv", package = "clinvaR")
    print(path)

    ## [1] "/Library/Frameworks/R.framework/Versions/3.4/Resources/library/clinvaR/extdata/MacArthur_Gene_Lists/lists/lmm_hcm.tsv"

    gene_panel <- get_genes(path)
    print(gene_panel)

    ##  [1] "ACTC1"  "ACTN2"  "CSRP3"  "GLA"    "LAMP2"  "MYBPC3" "MYH7"  
    ##  [8] "MYL2"   "MYL3"   "MYOZ2"  "NEXN"   "PLN"    "PRKAG2" "PTPN11"
    ## [15] "RAF1"   "TNNC1"  "TNNI3"  "TNNT2"  "TPM1"   "TTR"

#### From clinvaR stored lists

The MacArthur Lab has assembled a compilation of gene lists
(<https://github.com/macarthur-lab/gene_lists>) that may be a useful
starting point. We have included all sets (and a few of our own) with
abbreviated descriptions. Please visit the GitHub repo for citation
information.

<table>
<colgroup>
<col width="19%" />
<col width="8%" />
<col width="8%" />
<col width="62%" />
</colgroup>
<thead>
<tr class="header">
<th>List</th>
<th>File Name</th>
<th>Count</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Universe</td>
<td>universe.tsv</td>
<td>18,991</td>
<td>Approved symbols for 18,991 protein-coding genes according to HGNC. All other lists are subsets.</td>
</tr>
<tr class="even">
<td>Hypertrophic cardiomyopathy</td>
<td>lmm_hcm.tsv</td>
<td>20</td>
<td>Gene testing panel for hypertrophic cardiomyopathy from the Laboratory of Molecular Medicine</td>
</tr>
<tr class="odd">
<td>ACMG secondary findings recommendations v2.0</td>
<td>acmg_59.tsv</td>
<td>59</td>
<td>The American College of Medical Genetics and Genomics has recommended that sequencing laboratories seek and report secondary findings in a minimum list of 59 genes.</td>
</tr>
<tr class="even">
<td>BROCA - Cancer Risk Panel</td>
<td>BROCA_Cancer_Risk_Panel.tsv</td>
<td>66</td>
<td>Suspected hereditary cancer predisposition, with a focus on breast or ovarian cancers. May co-occur with other cancer types (such as colorectal, endometrial, pancreatic, endocrine, or melanoma).</td>
</tr>
<tr class="odd">
<td>DNA Repair Genes, KangJ</td>
<td>DRG_KangJ.tsv</td>
<td>151</td>
<td>151 DNA repair genes from DNA repair pathways: ATM, BER, FA/HR, MMR, NHEJ, NER, TLS, XLR, RECQ, and other.</td>
</tr>
<tr class="even">
<td>ClinGen haploinsufficient genes</td>
<td>clingen_level3_genes_2015_02_27.tsv</td>
<td>221</td>
<td>Genes with sufficient evidence for dosage pathogenicity (level 3) as determined by the ClinGen Dosage Sensitivity Map</td>
</tr>
<tr class="odd">
<td>Genes with any disease association reported in ClinVar</td>
<td>clinvar_path_likelypath.tsv</td>
<td>3078</td>
<td>All gene symbols for which there is at least one variant with an assertion of pathogenic or likely pathogenic in ClinVar.</td>
</tr>
<tr class="even">
<td>FDA-approved drug targets</td>
<td>fda_approved_drug_targets.tsv</td>
<td>286</td>
<td>Genes whose protein products are known to be the mechanistic targets of FDA-approved drugs.</td>
</tr>
<tr class="odd">
<td>Drug targets by Nelson et al 2012</td>
<td>drug_targets_nelson.tsv</td>
<td>201</td>
<td>Drug targets according to Nelson et al 2012.</td>
</tr>
<tr class="even">
<td>X-linked ClinVar genes</td>
<td>x-linked_clinvar.tsv</td>
<td>61</td>
<td>X chromosome genes in the August 6, 2015 ClinVar release that have at least 3 reportedly pathogenic, non-conflicted variants in ClinVar with at least one submitter other than OMIM or GeneReviews.</td>
</tr>
<tr class="odd">
<td>All dominant genes</td>
<td>all_ad.tsv</td>
<td>709</td>
<td>Currently the union of the Berg and Blekhman autosomal dominant OMIM gene lists, may add more lists later.</td>
</tr>
<tr class="even">
<td>All recessive genes</td>
<td>all_ar.tsv</td>
<td>1183</td>
<td>Currently the union of the Berg and Blekhman autosomal recessive OMIM gene lists, may add more lists later.</td>
</tr>
<tr class="odd">
<td>Essential in culture</td>
<td>core_essentials_hart.tsv</td>
<td>285</td>
<td>Genes deemed essential in multiple cultured cell lines based on shRNA screen data</td>
</tr>
<tr class="even">
<td>Essential in mice</td>
<td>mgi_essential.tsv</td>
<td>2,454</td>
<td>Genes where homozygous knockout in mice results in pre-, peri- or post-natal lethality.</td>
</tr>
<tr class="odd">
<td>Genes nearest to GWAS peaks</td>
<td>gwascatalog.tsv</td>
<td>3,762</td>
<td>Closest gene 3' and 5' of GWAS hits in the NHGRI GWAS catalog as of Feb 9, 2015</td>
</tr>
<tr class="even">
<td>Kinases</td>
<td>kinases.tsv</td>
<td>351</td>
<td>From UniProt's <a href="http://www.uniprot.org/docs/pkinfam">pkinfam list</a></td>
</tr>
<tr class="odd">
<td>G-protein-coupled receptors</td>
<td>gpcr.tsv</td>
<td>1705</td>
<td>GPCR list from <a href="http://www.guidetopharmacology.org/GRAC/GPCRListForward?class=A">guidetopharmacology.org</a></td>
</tr>
<tr class="even">
<td>Natural product targets</td>
<td>natural_product_targets.tsv</td>
<td>37</td>
<td>List of hand-curated targets of natural products</td>
</tr>
</tbody>
</table>

These are easily imported by directly referring to the file name.

    gene_panel <- get_genes('lmm_hcm.tsv')
    print(gene_panel)

    ##  [1] "ACTC1"  "ACTN2"  "CSRP3"  "GLA"    "LAMP2"  "MYBPC3" "MYH7"  
    ##  [8] "MYL2"   "MYL3"   "MYOZ2"  "NEXN"   "PLN"    "PRKAG2" "PTPN11"
    ## [15] "RAF1"   "TNNC1"  "TNNI3"  "TNNT2"  "TPM1"   "TTR"

The last way is to search OMIM for related genes.

    #omim_table <- read.table('clinvaR/storage/gene_lists/other_data/omim.full.tsv', 
    #                         sep = '\t', fill = T, header = T)
    #search_terms <- c("Cardiomyopathy")

Then, you can download variants of interest from 1000 Genomes

    #download_output <- download_1000g(genes = gene_panel)
    #final_vcf <- import_file_1000g(genes = gene_panel)

Import ClinVar VCF
------------------

This must be merged with some version of ClinVar. This can be done
either using the latest saved version (July 5, 2017), or by specifying a
date from the list.

    get_date_list()

    ##  [1] "2012-06-16" "2012-10-26" "2012-11-05" "2012-12-31" "2013-01-14"
    ##  [6] "2013-01-18" "2013-02-26" "2013-05-06" "2013-08-08" "2013-09-05"
    ## [11] "2013-09-30" "2013-11-05" "2013-12-03" "2013-12-30" "2014-02-11"
    ## [16] "2014-03-03" "2014-04-01" "2014-04-30" "2014-06-04" "2014-07-02"
    ## [21] "2014-08-07" "2014-09-02" "2014-09-29" "2014-11-05" "2014-12-02"
    ## [26] "2015-01-06" "2015-02-03" "2015-03-05" "2015-03-30" "2015-05-04"
    ## [31] "2015-06-03" "2015-06-29" "2015-08-04" "2015-09-01" "2015-09-29"
    ## [36] "2015-11-02" "2015-12-01" "2016-01-04" "2016-02-03" "2016-03-02"
    ## [41] "2016-04-05" "2016-05-02" "2016-05-31" "2016-07-05" "2016-08-02"
    ## [46] "2016-08-31" "2016-10-03" "2016-11-01" "2016-11-28" "2017-01-04"
    ## [51] "2017-01-30" "2017-02-28" "2017-04-04" "2017-05-01" "2017-05-30"
    ## [56] "2017-07-05"

Use `get_clinvar()` to collect from a certain date (if date is not
present, the closest date is used).

    clinvar <- get_clinvar("2017-07-21")
    clinvar[1:4,] #%>% select(-VAR_ID, -CLNDBN)

    ##         VAR_ID CHROM    POS REF ALT          ID CLNSIG
    ## 1 1_949523_C_T     1 949523   C   T rs786201005      5
    ## 2 1_949608_G_A     1 949608   G   A      rs1921      2
    ## 3 1_949739_G_T     1 949739   G   T rs672601312      5
    ## 4 1_955563_G_C     1 955563   G   C rs539283387      3
    ##                                                 CLNDBN
    ## 1 Immunodeficiency_38_with_basal_ganglia_calcification
    ## 2                                        not_specified
    ## 3 Immunodeficiency_38_with_basal_ganglia_calcification
    ## 4                                        not_specified
    ##   pathogenic_incl_conflicts pathogenic_no_conflicts
    ## 1                      TRUE                    TRUE
    ## 2                     FALSE                   FALSE
    ## 3                      TRUE                    TRUE
    ## 4                     FALSE                   FALSE

    clinvar <- get_clinvar()
    clinvar[1:4,] #%>% select(-VAR_ID, -CLNDBN)

    ##         VAR_ID CHROM    POS REF ALT          ID CLNSIG
    ## 1 1_949523_C_T     1 949523   C   T rs786201005      5
    ## 2 1_949608_G_A     1 949608   G   A      rs1921      2
    ## 3 1_949739_G_T     1 949739   G   T rs672601312      5
    ## 4 1_955563_G_C     1 955563   G   C rs539283387      3
    ##                                                 CLNDBN
    ## 1 Immunodeficiency_38_with_basal_ganglia_calcification
    ## 2                                        not_specified
    ## 3 Immunodeficiency_38_with_basal_ganglia_calcification
    ## 4                                        not_specified
    ##   pathogenic_incl_conflicts pathogenic_no_conflicts
    ## 1                      TRUE                    TRUE
    ## 2                     FALSE                   FALSE
    ## 3                      TRUE                    TRUE
    ## 4                     FALSE                   FALSE

Then, we merge clinvar with the genes of interest.

    #merged_1000g <- merge_clinvar_1000g(clinvar = get_clinvar(Sys.Date()), 
    #                    vcf = import_file_1000g(gene_panel))
    #merged_1000g[1:5,1:5]

We can also slice clinvar itself to collect only relevant entries.

    #clinvar$CHROM %>% paste0(str_pad(clinvar$POS, 10, side = 'left', pad = "0")) %>% as.numeric %>% 
    #  between(start, stop) -> filter_condition
    #clinvar <- filter(clinvar, filter_condition)
