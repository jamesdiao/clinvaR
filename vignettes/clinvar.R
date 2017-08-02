## ----setup, include = F--------------------------------------------------
knitr::opts_knit$set(root.dir = ".");
knitr::opts_chunk$set(echo = T, eval = T, cache = T, warning = F, message = F)

## ---- cache = T----------------------------------------------------------
hcm_vcf <- annotate_1000g(vcf = vcf, clinvar = newest_clinvar, conflicts = TRUE)
hcm_vcf[1:3,1:12]

## ---- fig.width = 6.5, fig.height = 4.5----------------------------------
plot(hcm_vcf, fraction = TRUE)
plot(hcm_vcf, fraction = FALSE)

## ---- fig.width = 6.5, fig.height = 4.5----------------------------------
plot(vcf, fraction = FALSE)

