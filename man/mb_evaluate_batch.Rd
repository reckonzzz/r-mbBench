\name{mb_evaluate_batch}
\alias{mb_evaluate_batch}
\title{Evaluate batch-correction results}
\usage{
mb_evaluate_batch(
  benchmark,
  input,
  output_dir = NULL,
  overwrite = FALSE,
  prevalence_min = 0.1,
  maximum_relative_abundance_min = 0.002,
  calculate_scores = FALSE
)
}
\arguments{
\item{benchmark}{An \code{mbbench_dataset}, loaded \code{mbbench_benchmark}, or bundle directory.}
\item{input}{A corrected sample-by-taxon data frame, matrix, CSV path, or named list of those values.}
\item{output_dir}{Optional directory for evaluation tables.}
\item{overwrite}{Whether to replace a nonempty output directory.}
\item{prevalence_min}{Observed prevalence threshold for structure metrics and optional scores.}
\item{maximum_relative_abundance_min}{Observed maximum-abundance threshold for structure metrics and optional scores.}
\item{calculate_scores}{Whether to calculate the separate Figure 4b-compatible component, category, and overall score table. This requires the R package \code{lisi} and the mbBench Python backend.}
}
\value{An \code{mbbench_batch_evaluation} list with method- and sample-level metrics, optional scores, applicability, resolved configuration, output directory, and files.}
\description{
Aligns each corrected matrix to the benchmark by sample and taxon labels. Paired
recovery metrics use every taxon; feature thresholds apply to structure metrics
and to LISI and classifier inputs when optional scores are requested.
}
