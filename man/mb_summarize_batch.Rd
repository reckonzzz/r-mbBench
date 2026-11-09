\name{mb_summarize_batch}
\alias{mb_summarize_batch}
\title{Summarize repeated batch evaluations}
\usage{
mb_summarize_batch(input, output_dir = NULL, overwrite = FALSE)
}
\arguments{
\item{input}{An evaluation, metrics data frame or CSV, directory, or list.}
\item{output_dir}{Optional summary directory.}
\item{overwrite}{Whether to replace a nonempty output directory.}
}
\value{An \code{mbbench_batch_summary}. \code{replicate_metrics} and optional
\code{replicate_scores} are wide tables with one row per method, scenario, and
seed. \code{metric_summary} and optional \code{score_summary} are long tables
with one row per method, scenario, and metric, including the repetition count,
mean, and standard deviation. Score fields are \code{NULL} when no
\code{batch_scores.csv} files are present.}
\description{Summarizes numeric batch-correction endpoints across benchmark repetitions and keeps optional publication scores in a separate table.}
