\name{mb_summarize_da}
\alias{mb_summarize_da}
\title{Summarize repeated DA evaluations}
\usage{
mb_summarize_da(
  input,
  output_dir = NULL,
  overwrite = FALSE,
  bootstrap_repetitions = 2000L,
  seed = 42L
)
}
\arguments{
\item{input}{An evaluation, metrics data frame or CSV, directory, or list.}
\item{output_dir}{Optional summary directory.}
\item{overwrite}{Whether to replace a nonempty output directory.}
\item{bootstrap_repetitions}{Number of bootstrap repetitions.}
\item{seed}{Bootstrap seed.}
}
\value{One summary row per method and scenario.}
\description{Summarizes seed-level FDP, power, type-I and type-II error, and FWER.}
