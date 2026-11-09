\name{mb_da_result}
\alias{mb_da_result}
\title{Build a complete DA decision table}
\usage{
mb_da_result(benchmark, significant_taxa, tested_taxa = NULL)
}
\arguments{
\item{benchmark}{An mbBench benchmark object or bundle directory.}
\item{significant_taxa}{Character vector of taxa called significant.}
\item{tested_taxa}{Optional character vector of tested taxa; defaults to all taxa.}
}
\value{A data frame with one row per benchmark taxon and logical \code{tested} and \code{significant} columns.}
\description{Creates a label-aligned DA input without manually filling non-significant taxa.}
