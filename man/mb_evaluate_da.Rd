\name{mb_evaluate_da}
\alias{mb_evaluate_da}
\title{Evaluate differential-abundance decisions}
\usage{
mb_evaluate_da(benchmark, input, output_dir = NULL, overwrite = FALSE)
}
\arguments{
\item{benchmark}{An mbBench benchmark object or bundle directory.}
\item{input}{A data frame, named logical vector, CSV path, or named list of method results.}
\item{output_dir}{Optional directory for evaluation tables.}
\item{overwrite}{Whether to replace a nonempty output directory.}
}
\value{An \code{mbbench_da_evaluation} list with metrics, taxon-level results, applicability, resolved configuration, output directory, and files.}
\description{
Evaluates one Boolean significance decision per benchmark taxon using the
bundle's declared truth applicability.
}
