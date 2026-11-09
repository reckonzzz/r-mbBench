\name{mb_evaluate_da_many}
\alias{mb_evaluate_da_many}
\title{Evaluate many DA benchmark repetitions}
\usage{
mb_evaluate_da_many(
  manifest,
  output_dir = NULL,
  overwrite = FALSE,
  bootstrap_repetitions = 2000L,
  seed = 42L
)
}
\arguments{
\item{manifest}{A data frame or CSV with \code{benchmark}, \code{input}, and \code{method} columns.}
\item{output_dir}{Optional evaluation directory.}
\item{overwrite}{Whether to replace a nonempty output directory.}
\item{bootstrap_repetitions}{Number of summary bootstrap repetitions.}
\item{seed}{Bootstrap seed.}
}
\value{An \code{mbbench_da_evaluation} containing repetition metrics, taxon results, applicability, and a summary.}
\description{Evaluates every manifest row and reports empirical FDR as the mean repetition-level FDP.}
