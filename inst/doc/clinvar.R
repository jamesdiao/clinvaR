## ----setup, include = F--------------------------------------------------
knitr::opts_knit$set(root.dir = ".");
knitr::opts_chunk$set(echo = T, eval = T, cache = T, warning = F, message = F)

## ---- fig.width = 6.5, fig.height = 4.5----------------------------------
over_time(vcf, verbose = FALSE) -> output
plot(output$Frequencies)
plot(output$Submissions)
plot(output$Significances)
#output$Frequencies %>% kable
#output$Submissions %>% kable

## ---- fig.width = 6.5, fig.height = 4.5----------------------------------
over_time(vcf, prevalence = 0.002, case_freq = 0.4, verbose = FALSE) -> output
plot(output$Penetrance)

