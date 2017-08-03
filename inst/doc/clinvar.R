## ----setup, include = F--------------------------------------------------
knitr::opts_knit$set(root.dir = ".");
knitr::opts_chunk$set(echo = T, eval = T, cache = T, warning = F, message = F)

## ---- fig.width = 6.5, fig.height = 4.5----------------------------------
over_time(vcf, verbose = FALSE) -> output
#output$Frequencies %>% kable
plot(output$Frequencies)
#output$Submissions %>% kable
plot(output$Submissions)

## ---- fig.width = 6.5, fig.height = 4.5----------------------------------
over_time(vcf, prevalence = 0.002, case_freq = 0.4, verbose = FALSE) -> output
#output$Frequencies %>% kable
plot(output$Penetrance)
#output$Submissions %>% kable

