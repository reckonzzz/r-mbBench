#' @export
print.mbbench_backend_status <- function(x, ...) {
  cat("<mbbench_backend_status>\n")
  cat(sprintf("  configured: %s\n", if (isTRUE(x$configured)) "yes" else "no"))
  cat(sprintf("  initialized: %s\n", if (isTRUE(x$initialized)) "yes" else "no"))
  available <- if (!isTRUE(x$checked)) "not checked" else if (isTRUE(x$available)) "yes" else "no"
  cat(sprintf("  available: %s\n", available))
  if (!is.null(x$python)) cat(sprintf("  python: %s\n", x$python))
  if (!is.null(x$version)) cat(sprintf("  version: %s\n", x$version))
  invisible(x)
}

#' @export
print.mbbench_dataset <- function(x, ...) {
  n_samples <- nrow(x$counts)
  n_taxa <- ncol(x$counts) - 1L
  cat(sprintf("<mbbench_dataset> %s samples x %s taxa\n", n_samples, n_taxa))
  cat(sprintf("  path: %s\n", x$path))
  cat(sprintf("  tables: %s\n", paste(names(x)[vapply(x, is.data.frame, logical(1))], collapse = ", ")))
  invisible(x)
}

#' @export
print.mbbench_model <- function(x, ...) {
  cat(sprintf("<mbbench_model> %s\n", x$model_type))
  cat(sprintf("  path: %s\n", x$path))
  cat(sprintf("  device: %s\n", x$device))
  cat(sprintf("  conditions: %s\n", paste(x$condition_fields, collapse = ", ")))
  invisible(x)
}

#' @export
print.mbbench_bank <- function(x, ...) {
  cat(sprintf("<mbbench_bank> %s samples x %s taxa\n", x$n_samples, x$n_taxa))
  cat(sprintf("  path: %s\n", x$path))
  cat(sprintf("  conditions: %s\n", paste(x$condition_fields, collapse = ", ")))
  invisible(x)
}

.mb_print_generated <- function(x, label) {
  cat(sprintf("<%s> estimated memory %.1f MiB\n", label, x$estimated_memory_mb))
  cat(sprintf("  loaded: %s\n", if (isTRUE(x$loaded)) "yes" else "no"))
  if (!is.null(x$output_dir)) cat(sprintf("  output: %s\n", x$output_dir))
  tables <- names(x)[vapply(x, is.data.frame, logical(1))]
  if (length(tables)) cat(sprintf("  in-memory tables: %s\n", paste(tables, collapse = ", ")))
  if (length(x$files)) cat(sprintf("  files: %s\n", paste(names(x$files), collapse = ", ")))
  invisible(x)
}

#' @export
print.mbbench_generation <- function(x, ...) {
  .mb_print_generated(x, "mbbench_generation")
}

#' @export
print.mbbench_benchmark <- function(x, ...) {
  .mb_print_generated(x, "mbbench_benchmark")
}

#' @export
print.mbbench_batch_evaluation <- function(x, ...) {
  cat(sprintf("<mbbench_batch_evaluation> %s method(s)\n", nrow(x$metrics)))
  if (!is.null(x$scores)) cat(sprintf("  score rows: %s\n", nrow(x$scores)))
  if (!is.null(x$output_dir)) cat(sprintf("  output: %s\n", x$output_dir))
  invisible(x)
}

#' @export
print.mbbench_batch_summary <- function(x, ...) {
  cat(sprintf("<mbbench_batch_summary> %s metric rows\n", nrow(x$metric_summary)))
  if (!is.null(x$score_summary)) {
    cat(sprintf("  score rows: %s\n", nrow(x$score_summary)))
  }
  if (!is.null(x$output_dir)) cat(sprintf("  output: %s\n", x$output_dir))
  invisible(x)
}

#' @export
print.mbbench_da_evaluation <- function(x, ...) {
  cat(sprintf("<mbbench_da_evaluation> %s run(s)\n", nrow(x$metrics)))
  if (!is.null(x$summary)) cat(sprintf("  summary rows: %s\n", nrow(x$summary)))
  if (!is.null(x$output_dir)) cat(sprintf("  output: %s\n", x$output_dir))
  invisible(x)
}
