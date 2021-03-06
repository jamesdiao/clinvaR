% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_clinvar.R
\name{get_closest_date}
\alias{get_closest_date}
\title{Get Closest ClinVar Version Dates}
\usage{
get_date_list() 
get_date_list(dates)
}
\arguments{
\item{file}{character; input date in the format \%y-\%m-\%d. Ex: 2017-07-05.}
}
\description{
This function returns the closest ClinVar date to the input date. 
If there is no input, it defaults to the current date.
}
