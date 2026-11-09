.mb_eval_bundle <- function(benchmark) {
  if (is.character(benchmark) && length(benchmark) == 1L) {
    return(mb_dataset(path = benchmark))
  }
  if (!is.list(benchmark)) {
    .mb_abort("benchmark must be an mbBench benchmark object or bundle directory.")
  }
  required <- c("counts", "relative_abundance", "metadata", "feature_truth",
                "batch_neutral_reference", "config")
  missing <- required[vapply(required, function(name) is.null(benchmark[[name]]), logical(1))]
  if (length(missing)) {
    .mb_abort("benchmark is not loaded; missing: %s.", paste(missing, collapse = ", "))
  }
  benchmark
}

.mb_eval_identity <- function(bundle) {
  scenario <- bundle$config$preset %||% bundle$config$scenario_id %||% "custom"
  seed <- bundle$config$runtime$seed %||% NA_integer_
  list(scenario = as.character(scenario), seed = as.integer(seed))
}

.mb_eval_bundle_path <- function(bundle) {
  path <- bundle$path %||% bundle$output_dir %||% NULL
  if (is.null(path)) NULL else as.character(path)
}

.mb_eval_frame <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    return(utils::read.csv(x, check.names = FALSE, stringsAsFactors = FALSE))
  }
  if (is.matrix(x)) return(as.data.frame(x, check.names = FALSE))
  if (is.data.frame(x)) return(x)
  .mb_abort("input must be a data.frame, matrix, or CSV path.")
}

.mb_close_rows <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!length(x) || any(!is.finite(x)) || any(x < 0)) {
    .mb_abort("Abundance input must contain finite nonnegative values.")
  }
  totals <- rowSums(x)
  if (any(totals <= 0)) .mb_abort("Every abundance row must have positive mass.")
  x / totals
}

.mb_clr <- function(x, pseudocount) {
  logged <- log(.mb_close_rows(x) + pseudocount)
  logged - rowMeans(logged)
}

.mb_method_inputs <- function(input) {
  if (is.list(input) && !is.data.frame(input)) {
    methods <- trimws(names(input) %||% character())
    if (!length(input) || length(methods) != length(input) ||
        any(!nzchar(methods)) || anyDuplicated(methods)) {
      .mb_abort("A multi-method input must be a uniquely named list.")
    }
    names(input) <- methods
    return(input)
  }
  list(method = input)
}

.mb_bind_rows <- function(rows) {
  fields <- unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned <- lapply(rows, function(row) {
    missing <- setdiff(fields, names(row))
    if (length(missing)) row[missing] <- NA
    row[fields]
  })
  result <- do.call(rbind, aligned)
  rownames(result) <- NULL
  result
}

.mb_batch_matrix <- function(input, sample_ids, taxa) {
  frame <- .mb_eval_frame(input)
  if ("sample_id" %in% names(frame)) {
    ids <- as.character(frame$sample_id)
    frame$sample_id <- NULL
  } else {
    ids <- rownames(frame)
    if (is.null(ids) || identical(ids, as.character(seq_len(nrow(frame))))) {
      .mb_abort("Batch input requires a sample_id column or informative row names.")
    }
  }
  if (anyDuplicated(ids)) .mb_abort("Batch input sample IDs must be unique.")
  missing_samples <- setdiff(sample_ids, ids)
  extra_samples <- setdiff(ids, sample_ids)
  missing_taxa <- setdiff(taxa, names(frame))
  extra_taxa <- setdiff(names(frame), taxa)
  if (length(missing_samples) || length(extra_samples)) {
    .mb_abort("Batch input sample IDs do not match the benchmark.")
  }
  if (length(missing_taxa) || length(extra_taxa)) {
    .mb_abort("Batch input taxa do not match the benchmark.")
  }
  rownames(frame) <- ids
  .mb_close_rows(frame[sample_ids, taxa, drop = FALSE])
}

.mb_paired_distances <- function(candidate, reference, pseudocount) {
  candidate <- .mb_close_rows(candidate)
  reference <- .mb_close_rows(reference)
  if (!identical(dim(candidate), dim(reference))) {
    .mb_abort("Candidate and reference matrices must have the same shape.")
  }
  error <- .mb_clr(candidate, pseudocount) - .mb_clr(reference, pseudocount)
  squared <- rowSums(error ^ 2)
  list(
    aitchison = sqrt(squared),
    clr_rms = sqrt(squared / ncol(candidate)),
    bray_curtis = rowSums(abs(candidate - reference)) / 2,
    hellinger = sqrt(rowSums((sqrt(candidate) - sqrt(reference)) ^ 2)) / sqrt(2)
  )
}

.mb_factor_matrix <- function(x) {
  text <- as.character(x)
  levels <- sort(unique(text))
  if (length(levels) < 2L) return(matrix(numeric(), nrow = length(text), ncol = 0L))
  vapply(levels[-1L], function(level) as.numeric(text == level), numeric(length(text)))
}

.mb_design <- function(...) {
  variables <- list(...)
  columns <- list(rep(1, length(variables[[1L]])))
  for (variable in variables) {
    encoded <- .mb_factor_matrix(variable)
    if (ncol(encoded)) {
      columns <- c(columns, lapply(seq_len(ncol(encoded)), function(index) encoded[, index]))
    }
  }
  do.call(cbind, columns)
}

.mb_projection <- function(design) {
  decomposition <- qr(design)
  if (!decomposition$rank) return(matrix(0, nrow(design), nrow(design)))
  basis <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
  basis %*% t(basis)
}

.mb_marginal_r2 <- function(distance, reduced, full) {
  n <- nrow(distance)
  center <- diag(n) - matrix(1 / n, n, n)
  gower <- -0.5 * center %*% (distance ^ 2) %*% center
  total <- sum(diag(gower))
  if (!is.finite(total) || total <= 0) return(NA_real_)
  effect <- sum(diag((.mb_projection(full) - .mb_projection(reduced)) %*% gower))
  max(0, effect) / total
}

.mb_stratified_asw <- function(embedding, batch, group) {
  values <- weights <- numeric()
  distance <- as.matrix(stats::dist(embedding))
  for (level in unique(as.character(group))) {
    keep <- which(as.character(group) == level)
    labels <- as.character(batch[keep])
    counts <- table(labels)
    if (length(counts) < 2L || min(counts) < 2L) {
      return(list(value = NA_real_, reason = sprintf(
        "batch ASW is not estimable in biological group '%s'", level
      )))
    }
    scores <- vapply(seq_along(keep), function(i) {
      own <- keep[i]
      same <- keep[labels == labels[i] & keep != own]
      other_labels <- setdiff(unique(labels), labels[i])
      a <- mean(distance[own, same])
      b <- min(vapply(other_labels, function(other) {
        mean(distance[own, keep[labels == other]])
      }, numeric(1)))
      if (max(a, b) == 0) 0 else (b - a) / max(a, b)
    }, numeric(1))
    values <- c(values, mean(scores))
    weights <- c(weights, length(keep))
  }
  list(value = stats::weighted.mean(values, weights), reason = NULL)
}

.mb_structure_context <- function(counts, corrected, reference, metadata,
                                  prevalence_min, maximum_relative_abundance_min) {
  if (!all(c("group", "batch") %in% names(metadata))) {
    .mb_abort("Benchmark metadata requires group and batch columns.")
  }
  observed <- .mb_close_rows(counts)
  active <- colMeans(counts > 0) >= prevalence_min &
    apply(observed, 2L, max) >= maximum_relative_abundance_min
  result <- list(
    structure_prevalence_min = prevalence_min,
    structure_maximum_relative_abundance_min = maximum_relative_abundance_min,
    structure_n_active_taxa = sum(active),
    structure_metrics_supported = FALSE,
    structure_metrics_reason = NA_character_
  )
  if (sum(active) < 2L) {
    result$structure_metrics_reason <- "fewer than two taxa pass the R2/ASW feature filter"
    return(list(metrics = result, context = NULL))
  }
  active_counts <- counts[, active, drop = FALSE]
  active_corrected <- corrected[, active, drop = FALSE]
  active_reference <- reference[, active, drop = FALSE]
  if (any(rowSums(active_counts) <= 0) || any(rowSums(active_corrected) <= 0) ||
      any(rowSums(active_reference) <= 0)) {
    result$structure_metrics_reason <-
      "the R2/ASW feature filter leaves at least one zero-mass sample"
    return(list(metrics = result, context = NULL))
  }
  pseudocount <- max(0.5 / stats::median(rowSums(active_counts)), 1e-12)
  stages <- list(
    observed = .mb_clr(active_counts, pseudocount),
    corrected = .mb_clr(active_corrected, pseudocount),
    batch_neutral_reference = .mb_clr(active_reference, pseudocount)
  )
  batch <- as.character(metadata$batch)
  group <- as.character(metadata$group)
  list(metrics = result, context = list(
    stages = stages, batch = batch, group = group, pseudocount = pseudocount
  ))
}

.mb_shared_embeddings <- function(stages) {
  keep <- apply(stages$observed, 2L, stats::sd) > 0
  components <- min(20L, nrow(stages$observed) - 1L, sum(keep))
  if (components < 2L) {
    return(list(
      embeddings = NULL,
      components = components,
      reason = "fewer than two nonconstant observed CLR features remain for ASW"
    ))
  }
  pca <- stats::prcomp(stages$observed[, keep, drop = FALSE], center = TRUE,
                       scale. = FALSE, rank. = components)
  list(
    embeddings = lapply(stages, function(stage) {
      stats::predict(pca, stage[, keep, drop = FALSE])[
        , seq_len(components), drop = FALSE
      ]
    }),
    components = components,
    reason = NULL
  )
}

.mb_structure_metrics <- function(counts, corrected, reference, metadata,
                                  prevalence_min, maximum_relative_abundance_min) {
  prepared <- .mb_structure_context(
    counts, corrected, reference, metadata,
    prevalence_min, maximum_relative_abundance_min
  )
  result <- prepared$metrics
  if (is.null(prepared$context)) return(result)
  stages <- prepared$context$stages
  batch <- prepared$context$batch
  group <- prepared$context$group
  biology_design <- .mb_design(group)
  batch_design <- .mb_design(batch)
  full_design <- .mb_design(group, batch)
  for (stage in names(stages)) {
    distance <- as.matrix(stats::dist(stages[[stage]]))
    result[[paste0(stage, "_batch_aitchison_marginal_r2")]] <-
      .mb_marginal_r2(distance, biology_design, full_design)
    result[[paste0(stage, "_biology_aitchison_marginal_r2")]] <-
      .mb_marginal_r2(distance, batch_design, full_design)
  }
  observed_biology <- result$observed_biology_aitchison_marginal_r2
  corrected_biology <- result$corrected_biology_aitchison_marginal_r2
  result$biology_aitchison_r2_change <- corrected_biology - observed_biology
  observed_batch <- result$observed_batch_aitchison_marginal_r2
  corrected_batch <- result$corrected_batch_aitchison_marginal_r2
  result$batch_aitchison_r2_change <- observed_batch - corrected_batch
  embedded <- .mb_shared_embeddings(stages)
  if (is.null(embedded$embeddings)) {
    result$structure_metrics_reason <- embedded$reason
    return(result)
  }
  for (stage in names(embedded$embeddings)) {
    embedding <- embedded$embeddings[[stage]]
    asw <- .mb_stratified_asw(embedding, batch, group)
    result[[paste0(stage, "_batch_asw_stratified")]] <- asw$value
    if (!is.null(asw$reason)) result$structure_metrics_reason <- asw$reason
  }
  result$batch_asw_change <- result$observed_batch_asw_stratified -
    result$corrected_batch_asw_stratified
  result$structure_relative_pseudocount <- prepared$context$pseudocount
  result$structure_pca_components <- embedded$components
  result$structure_metrics_supported <- is.na(result$structure_metrics_reason)
  result
}

.mb_score_mean <- function(values) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values)) mean(values) else NA_real_
}

.mb_bounded_similarity <- function(observed, corrected) {
  if (length(observed) != 1L || length(corrected) != 1L ||
      !is.finite(observed) || !is.finite(corrected)) return(NA_real_)
  scale <- max(abs(observed), abs(corrected), .Machine$double.eps)
  max(0, min(1, 1 - abs(corrected - observed) / scale))
}

.mb_removal_fraction <- function(observed, corrected) {
  if (length(observed) != 1L || length(corrected) != 1L ||
      !is.finite(observed) || !is.finite(corrected) ||
      observed <= .Machine$double.eps) return(NA_real_)
  max(0, min(1, (observed - corrected) / observed))
}

.mb_mixing_improvement <- function(observed, corrected) {
  if (length(observed) != 1L || length(corrected) != 1L ||
      !is.finite(observed) || !is.finite(corrected)) return(NA_real_)
  denominator <- 1 - observed
  if (denominator <= .Machine$double.eps) return(NA_real_)
  max(0, min(1, (corrected - observed) / denominator))
}

.mb_lisi_expected_adjust <- function(values, labels) {
  probabilities <- as.numeric(table(labels)) / length(labels)
  expected <- 1 / sum(probabilities ^ 2)
  denominator <- expected - 1
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
    return(rep(NA_real_, length(values)))
  }
  pmin(1, pmax(0, (values - 1) / denominator))
}

.mb_lisi_value <- function(values, labels, compute_lisi) {
  labels <- factor(labels)
  if (nlevels(labels) < 2L || min(table(labels)) < 2L || nrow(values) < 7L) {
    return(NA_real_)
  }
  perplexity <- min(30, max(2, floor((nrow(values) - 1) / 3)))
  result <- tryCatch(
    compute_lisi(
      values,
      data.frame(label = labels),
      label_colnames = "label",
      perplexity = perplexity
    )$label,
    error = function(error) NULL
  )
  if (is.null(result)) return(NA_real_)
  adjusted <- .mb_lisi_expected_adjust(result, labels)
  adjusted <- adjusted[is.finite(adjusted)]
  if (length(adjusted)) stats::median(adjusted) else NA_real_
}

.mb_standard_lisi <- function(embeddings, metadata) {
  if (!requireNamespace("lisi", quietly = TRUE)) {
    .mb_abort(paste0(
      "calculate_scores = TRUE requires the R package 'lisi'. Install it with ",
      "remotes::install_github(\"immunogenomics/LISI\")."
    ))
  }
  compute_lisi <- getExportedValue("lisi", "compute_lisi")
  stage_metric <- function(values) {
    strata <- as.character(metadata$group)
    levels <- unique(strata)
    batch_lisi <- vapply(levels, function(level) {
      keep <- strata == level
      .mb_lisi_value(values[keep, , drop = FALSE], metadata$batch[keep], compute_lisi)
    }, numeric(1))
    weights <- vapply(levels, function(level) sum(strata == level), numeric(1))
    keep <- is.finite(batch_lisi)
    batch_value <- if (any(keep)) {
      stats::weighted.mean(batch_lisi[keep], weights[keep])
    } else {
      NA_real_
    }
    biology_lisi <- .mb_lisi_value(values, metadata$group, compute_lisi)
    list(
      batch_ilisi = batch_value,
      biology_separation = if (is.finite(biology_lisi)) 1 - biology_lisi else NA_real_
    )
  }
  list(
    stages = lapply(embeddings, stage_metric),
    versions = list(
      r_version = as.character(getRversion()),
      lisi_version = as.character(utils::packageVersion("lisi"))
    )
  )
}

.mb_cross_batch_classifier <- function(module, representation, metadata, seed) {
  result <- module$cross_batch_classifier_metrics(
    representation,
    as.character(metadata$group),
    as.character(metadata$batch),
    seed = as.integer(seed)
  )
  if (inherits(result, "python.builtin.dict")) result <- reticulate::py_to_r(result)
  result
}

.mb_classifier_preservation <- function(observed, corrected) {
  components <- vapply(c("auroc", "balanced_accuracy", "macro_f1"), function(metric) {
    .mb_bounded_similarity(observed[[metric]], corrected[[metric]])
  }, numeric(1))
  if (all(is.finite(components))) .mb_score_mean(components) else NA_real_
}

.mb_batch_score_row <- function(counts, corrected, reference, metadata, raw,
                                method, scenario, seed, prevalence_min,
                                maximum_relative_abundance_min,
                                classifier_module) {
  prepared <- .mb_structure_context(
    counts, corrected, reference, metadata,
    prevalence_min, maximum_relative_abundance_min
  )
  batch_applicable <- isTRUE(raw$batch_recovery_applicable)
  biology_applicable <- isTRUE(raw$biological_preservation_applicable)
  classifier_seed <- as.integer((if (is.na(seed)) 42L else seed) + 17L)
  empty_scores <- list(
    batch_global_score = NA_real_, batch_local_score = NA_real_,
    batch_effect_removal_score = NA_real_,
    biology_global_preservation_score = NA_real_,
    biology_local_preservation_score = NA_real_,
    biology_classifier_preservation_score = NA_real_,
    biological_signal_preservation_score = NA_real_, overall_score = NA_real_
  )
  base <- list(
    method = method, scenario = scenario, seed = seed,
    batch_recovery_applicable = batch_applicable,
    biological_preservation_applicable = biology_applicable,
    score_definition = "figure4b_publication_v1",
    classifier_seed = classifier_seed
  )
  if (is.null(prepared$context)) {
    return(list(
      row = c(base, list(score_supported = FALSE,
                         score_reason = prepared$metrics$structure_metrics_reason),
              empty_scores),
      versions = list()
    ))
  }
  embedded <- tryCatch(
    classifier_module$publication_lisi_embeddings(
      prepared$context$stages,
      seed = as.integer(if (is.na(seed)) 42L else seed)
    ),
    error = identity
  )
  if (inherits(embedded, "error")) {
    return(list(
      row = c(base, list(score_supported = FALSE,
                         score_reason = conditionMessage(embedded)),
              empty_scores),
      versions = list()
    ))
  }
  if (inherits(embedded, "python.builtin.dict")) embedded <- reticulate::py_to_r(embedded)

  lisi <- .mb_standard_lisi(embedded, metadata)
  classifiers <- list()
  if (biology_applicable) {
    classifiers$observed <- .mb_cross_batch_classifier(
      classifier_module, prepared$context$stages$observed, metadata, classifier_seed
    )
    classifiers$corrected <- .mb_cross_batch_classifier(
      classifier_module, prepared$context$stages$corrected, metadata, classifier_seed
    )
  }

  scores <- empty_scores
  if (batch_applicable) {
    scores$batch_global_score <- .mb_removal_fraction(
      raw$observed_batch_aitchison_marginal_r2,
      raw$corrected_batch_aitchison_marginal_r2
    )
    asw_observed <- 1 - abs(raw$observed_batch_asw_stratified)
    asw_corrected <- 1 - abs(raw$corrected_batch_asw_stratified)
    asw_score <- .mb_mixing_improvement(asw_observed, asw_corrected)
    batch_lisi_values <- vapply(
      lisi$stages, `[[`, numeric(1), "batch_ilisi"
    )
    lisi_score <- if (all(is.finite(batch_lisi_values))) {
      .mb_mixing_improvement(batch_lisi_values[["observed"]],
                             batch_lisi_values[["corrected"]])
    } else {
      NA_real_
    }
    scores$batch_local_score <- if (is.finite(asw_score) && is.finite(lisi_score)) {
      .mb_score_mean(c(asw_score, lisi_score))
    } else {
      NA_real_
    }
    scores$batch_effect_removal_score <- if (
      is.finite(scores$batch_global_score) && is.finite(scores$batch_local_score)
    ) {
      .mb_score_mean(c(scores$batch_global_score, scores$batch_local_score))
    } else {
      NA_real_
    }
  }
  if (biology_applicable) {
    scores$biology_global_preservation_score <- .mb_bounded_similarity(
      raw$observed_biology_aitchison_marginal_r2,
      raw$corrected_biology_aitchison_marginal_r2
    )
    biology_lisi_values <- vapply(
      lisi$stages, `[[`, numeric(1), "biology_separation"
    )
    scores$biology_local_preservation_score <- if (all(is.finite(biology_lisi_values))) {
      .mb_bounded_similarity(biology_lisi_values[["observed"]],
                             biology_lisi_values[["corrected"]])
    } else {
      NA_real_
    }
    scores$biology_classifier_preservation_score <- .mb_classifier_preservation(
      classifiers$observed, classifiers$corrected
    )
    biology_components <- c(
      scores$biology_global_preservation_score,
      scores$biology_local_preservation_score,
      scores$biology_classifier_preservation_score
    )
    scores$biological_signal_preservation_score <- if (all(is.finite(biology_components))) {
      .mb_score_mean(biology_components)
    } else {
      NA_real_
    }
  }
  if (is.finite(scores$batch_effect_removal_score) &&
      is.finite(scores$biological_signal_preservation_score)) {
    scores$overall_score <- 0.4 * scores$batch_effect_removal_score +
      0.6 * scores$biological_signal_preservation_score
  }
  expected <- character()
  if (batch_applicable) {
    expected <- c(expected, "batch_global_score", "batch_local_score",
                  "batch_effect_removal_score")
  }
  if (biology_applicable) {
    expected <- c(expected, "biology_global_preservation_score",
                  "biology_local_preservation_score",
                  "biology_classifier_preservation_score",
                  "biological_signal_preservation_score")
  }
  if (batch_applicable && biology_applicable) expected <- c(expected, "overall_score")
  missing <- expected[!vapply(scores[expected], is.finite, logical(1))]
  reasons <- character()
  if (length(classifiers)) {
    classifier_reasons <- vapply(classifiers, function(value) {
      value$reason %||% ""
    }, character(1))
    reasons <- classifier_reasons[nzchar(classifier_reasons)]
  }
  if (length(missing)) {
    reasons <- c(reasons, paste("unavailable score components:", paste(missing, collapse = ", ")))
  }
  list(
    row = c(base, list(score_supported = !length(missing),
                       score_reason = if (length(reasons)) paste(reasons, collapse = "; ") else NA_character_),
            scores),
    versions = lisi$versions
  )
}

.mb_eval_write <- function(output_dir, overwrite, files) {
  if (is.null(output_dir)) return(list(output_dir = NULL, files = list()))
  overwrite <- .mb_scalar_logical(overwrite, "overwrite")
  output <- .mb_normalize_path(output_dir, must_work = FALSE)
  if (file.exists(output) && !dir.exists(output) && !overwrite) {
    .mb_abort("output_dir is an existing file; use overwrite = TRUE to replace it.")
  }
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE)) &&
      !overwrite) {
    .mb_abort("output_dir is not empty; use overwrite = TRUE to replace it.")
  }
  parent <- dirname(output)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  staging <- tempfile(paste0(".", basename(output), ".staging-"), tmpdir = parent)
  dir.create(staging, recursive = TRUE, showWarnings = FALSE)
  on.exit({
    if (dir.exists(staging)) unlink(staging, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  written <- list()
  for (name in names(files)) {
    value <- files[[name]]
    path <- file.path(staging, name)
    if (grepl("\\.ya?ml$", name)) {
      yaml::write_yaml(value, path)
    } else if (grepl("\\.gz$", name)) {
      connection <- gzfile(path, "wt")
      tryCatch(utils::write.csv(value, connection, row.names = FALSE),
               finally = close(connection))
    } else {
      utils::write.csv(value, path, row.names = FALSE)
    }
    written[[sub("\\..*$", "", name)]] <- file.path(output, name)
  }
  if (file.exists(output)) {
    if (dir.exists(output)) unlink(output, recursive = TRUE, force = TRUE)
    else unlink(output, force = TRUE)
  }
  if (!file.rename(staging, output)) {
    .mb_abort("Could not install completed evaluation output in %s.", output)
  }
  list(output_dir = output, files = written)
}

#' Evaluate batch-correction results
#'
#' `input` is one corrected sample-by-taxon abundance matrix, or a named list
#' of matrices for several methods. Rows and taxa are aligned by labels.
#' Prevalence and maximum-abundance filters apply to structure metrics and the
#' optional score inputs; paired recovery metrics always use every benchmark taxon.
#'
#' @param benchmark An `mbbench_dataset`, loaded `mbbench_benchmark`, or bundle directory.
#' @param input A data frame, matrix, CSV path, or named list of those values.
#' @param output_dir Optional directory for evaluation tables.
#' @param overwrite Replace a nonempty output directory.
#' @param prevalence_min Observed prevalence threshold for structure metrics and optional scores.
#' @param maximum_relative_abundance_min Observed maximum-abundance threshold for structure metrics and optional scores.
#' @param calculate_scores Also calculate the separate publication score table.
#' @return An `mbbench_batch_evaluation` list containing `metrics`,
#'   `sample_metrics`, optional `scores`, `config`, `output_dir`, and `files`.
#' @export
mb_evaluate_batch <- function(benchmark, input, output_dir = NULL,
                              overwrite = FALSE, prevalence_min = 0.10,
                              maximum_relative_abundance_min = 0.002,
                              calculate_scores = FALSE) {
  bundle <- .mb_eval_bundle(benchmark)
  if (!is.logical(calculate_scores) || length(calculate_scores) != 1L ||
      is.na(calculate_scores)) {
    .mb_abort("calculate_scores must be TRUE or FALSE.")
  }
  prevalence_min <- .mb_scalar_number(prevalence_min, "prevalence_min", minimum = 0)
  if (prevalence_min > 1) .mb_abort("prevalence_min must not exceed 1.")
  maximum_relative_abundance_min <- .mb_scalar_number(
    maximum_relative_abundance_min, "maximum_relative_abundance_min", minimum = 0
  )
  taxa <- setdiff(names(bundle$counts), "sample_id")
  sample_ids <- as.character(bundle$counts$sample_id)
  counts <- as.matrix(bundle$counts[taxa])
  observed <- as.matrix(bundle$relative_abundance[taxa])
  reference <- as.matrix(bundle$batch_neutral_reference[taxa])
  pseudocount <- max(0.5 / stats::median(rowSums(counts)), 1e-12)
  identity <- .mb_eval_identity(bundle)
  method_inputs <- .mb_method_inputs(input)
  metrics <- sample_metrics <- scores <- list()
  score_versions <- list()
  classifier_module <- NULL
  if (isTRUE(calculate_scores)) {
    if (!requireNamespace("lisi", quietly = TRUE)) {
      .mb_abort(paste0(
        "calculate_scores = TRUE requires the R package 'lisi'. Install it with ",
        "remotes::install_github(\"immunogenomics/LISI\")."
      ))
    }
    .mb_ensure_backend()
    classifier_module <- reticulate::import("mbbench.batch_scores", convert = TRUE)
  }
  for (method in names(method_inputs)) {
    corrected <- .mb_batch_matrix(method_inputs[[method]], sample_ids, taxa)
    before <- .mb_paired_distances(observed, reference, pseudocount)
    after <- .mb_paired_distances(corrected, reference, pseudocount)
    unchanged <- abs(before$aitchison - after$aitchison) <=
      (1e-12 + 1e-10 * abs(after$aitchison))
    closer <- after$aitchison < before$aitchison & !unchanged
    farther <- after$aitchison > before$aitchison & !unchanged
    observed_rms <- sqrt(mean(before$aitchison ^ 2))
    corrected_rms <- sqrt(mean(after$aitchison ^ 2))
    pbrs_raw <- if (observed_rms <= .Machine$double.eps) NA_real_ else
      1 - corrected_rms ^ 2 / observed_rms ^ 2
    row <- list(
      method = method, scenario = identity$scenario, seed = identity$seed,
      batch_recovery_applicable = isTRUE(
        bundle$config$truth_applicability$batch_recovery
      ),
      biological_preservation_applicable = isTRUE(
        bundle$config$truth_applicability$biological_preservation
      ),
      relative_pseudocount = pseudocount,
      observed_rms_aitchison = observed_rms,
      corrected_rms_aitchison = corrected_rms,
      mean_paired_aitchison_improvement = mean(before$aitchison - after$aitchison),
      fraction_samples_closer = mean(closer),
      fraction_samples_unchanged = mean(unchanged),
      fraction_samples_farther = mean(farther),
      pbrs_raw = pbrs_raw,
      pbrs = if (is.na(pbrs_raw)) NA_real_ else max(0, min(1, pbrs_raw))
    )
    for (metric in c("clr_rms", "bray_curtis", "hellinger")) {
      row[[paste0("observed_mean_", metric)]] <- mean(before[[metric]])
      row[[paste0("corrected_mean_", metric)]] <- mean(after[[metric]])
    }
    row <- c(row, .mb_structure_metrics(
      counts, corrected, reference, bundle$metadata,
      prevalence_min, maximum_relative_abundance_min
    ))
    metrics[[method]] <- as.data.frame(row, check.names = FALSE)
    if (isTRUE(calculate_scores)) {
      scored <- .mb_batch_score_row(
        counts, corrected, reference, bundle$metadata, row,
        method, identity$scenario, identity$seed,
        prevalence_min, maximum_relative_abundance_min, classifier_module
      )
      scores[[method]] <- as.data.frame(scored$row, check.names = FALSE)
      if (length(scored$versions)) {
        if (length(score_versions) && !identical(score_versions, scored$versions)) {
          .mb_abort("R or LISI version changed during evaluation.")
        }
        score_versions <- scored$versions
      }
    }
    sample_metrics[[method]] <- data.frame(
      sample_id = sample_ids, method = method, scenario = identity$scenario,
      seed = identity$seed, observed_aitchison = before$aitchison,
      corrected_aitchison = after$aitchison,
      aitchison_improvement = before$aitchison - after$aitchison,
      observed_clr_rms = before$clr_rms, corrected_clr_rms = after$clr_rms,
      observed_bray_curtis = before$bray_curtis,
      corrected_bray_curtis = after$bray_curtis,
      observed_hellinger = before$hellinger,
      corrected_hellinger = after$hellinger,
      status = ifelse(closer, "closer", ifelse(farther, "farther", "unchanged")),
      check.names = FALSE
    )
  }
  metrics <- .mb_bind_rows(metrics)
  sample_metrics <- .mb_bind_rows(sample_metrics)
  scores <- if (isTRUE(calculate_scores)) .mb_bind_rows(scores) else NULL
  config <- list(
    schema_version = "mbbench_batch_evaluation_v1",
    benchmark = list(path = .mb_eval_bundle_path(bundle),
                     scenario = identity$scenario, seed = identity$seed),
    methods = names(method_inputs),
    paired_metrics_use_full_taxon_roster = TRUE,
    structure_feature_filter = list(
      prevalence_min = prevalence_min,
      maximum_relative_abundance_min = maximum_relative_abundance_min,
      applies_to = c(
        "aitchison_marginal_r2", "batch_asw",
        if (isTRUE(calculate_scores)) c("score_lisi", "score_classifier")
      ),
      does_not_apply_to = c("paired_recovery", "da_evaluation")
    ),
    truth_applicability = bundle$config$truth_applicability %||% list(),
    calculate_scores = isTRUE(calculate_scores)
  )
  if (isTRUE(calculate_scores)) {
    config$scores <- list(
      definition = "figure4b_publication_v1",
      file = "batch_scores.csv",
      r_version = score_versions$r_version %||% NULL,
      lisi_version = score_versions$lisi_version %||% NULL,
      classifier = "leave_one_batch_out_logistic_regression"
    )
  }
  output_files <- list(
    batch_metrics.csv = metrics,
    batch_sample_metrics.csv.gz = sample_metrics,
    evaluation.yaml = config
  )
  if (isTRUE(calculate_scores)) output_files$batch_scores.csv <- scores
  written <- .mb_eval_write(output_dir, overwrite, output_files)
  if (length(written$files)) {
    order <- c("batch_metrics", "batch_sample_metrics")
    labels <- c("metrics", "sample_metrics")
    if (isTRUE(calculate_scores)) {
      order <- c(order, "batch_scores")
      labels <- c(labels, "scores")
    }
    order <- c(order, "evaluation")
    labels <- c(labels, "config")
    written$files <- written$files[order]
    names(written$files) <- labels
  }
  structure(list(metrics = metrics, sample_metrics = sample_metrics, scores = scores,
                 applicability = bundle$config$truth_applicability %||% list(),
                 config = config, output_dir = written$output_dir, files = written$files),
            class = c("mbbench_batch_evaluation", "list"))
}

.mb_boolean <- function(x, name) {
  if (is.logical(x) && !anyNA(x)) return(x)
  text <- tolower(trimws(as.character(x)))
  valid <- text %in% c("true", "false", "1", "0")
  if (anyNA(text) || !all(valid)) .mb_abort("%s must contain only Boolean values.", name)
  text %in% c("true", "1")
}

.mb_da_frame <- function(input) {
  if (is.logical(input) && !is.null(names(input))) {
    return(data.frame(taxon = names(input), significant = unname(input)))
  }
  .mb_eval_frame(input)
}

.mb_align_da <- function(input, taxa, method) {
  frame <- .mb_da_frame(input)
  if (!all(c("taxon", "significant") %in% names(frame))) {
    .mb_abort("DA input requires taxon and significant columns.")
  }
  frame$taxon <- as.character(frame$taxon)
  if (anyDuplicated(frame$taxon)) .mb_abort("DA input taxa must be unique.")
  if (length(setdiff(frame$taxon, taxa)) || length(setdiff(taxa, frame$taxon))) {
    .mb_abort("DA input taxa do not match the benchmark.")
  }
  frame$significant <- .mb_boolean(frame$significant, "significant")
  frame$tested <- if ("tested" %in% names(frame)) .mb_boolean(frame$tested, "tested") else TRUE
  if (any(frame$significant & !frame$tested)) {
    .mb_abort("A taxon cannot be significant when tested is false.")
  }
  for (column in intersect(c("p_value", "q_value"), names(frame))) {
    original <- frame[[column]]
    suppressWarnings(frame[[column]] <- as.numeric(as.character(original)))
    if (any(!is.na(original) & is.na(frame[[column]]))) {
      .mb_abort("%s must contain numeric values or missing values.", column)
    }
    invalid <- !is.na(frame[[column]]) & (frame[[column]] < 0 | frame[[column]] > 1)
    if (any(invalid)) .mb_abort("%s must be missing or lie in [0, 1].", column)
  }
  if ("effect" %in% names(frame)) {
    original <- frame$effect
    suppressWarnings(frame$effect <- as.numeric(as.character(original)))
    if (any(!is.na(original) & is.na(frame$effect))) {
      .mb_abort("effect must contain numeric values or missing values.")
    }
    if (any(is.infinite(frame$effect), na.rm = TRUE)) {
      .mb_abort("effect must contain finite values or missing values.")
    }
  }
  frame <- merge(data.frame(taxon = taxa), frame, by = "taxon", sort = FALSE)
  frame <- frame[match(taxa, frame$taxon), , drop = FALSE]
  frame$method <- method
  frame[c("method", setdiff(names(frame), "method"))]
}

.mb_da_inputs <- function(input) {
  if (is.list(input) && !is.data.frame(input)) return(.mb_method_inputs(input))
  list(method = input)
}

.mb_da_truth <- function(feature_truth) {
  n <- nrow(feature_truth)
  abundance <- if ("biological_abundance_coefficient" %in% names(feature_truth))
    as.numeric(feature_truth$biological_abundance_coefficient) else rep(0, n)
  support <- if ("biological_support_coefficient" %in% names(feature_truth))
    as.numeric(feature_truth$biological_support_coefficient) else rep(0, n)
  direct <- if ("is_direct_biological_target" %in% names(feature_truth)) {
    .mb_boolean(feature_truth$is_direct_biological_target, "is_direct_biological_target")
  } else if ("biological_target" %in% names(feature_truth)) {
    .mb_boolean(feature_truth$biological_target, "biological_target")
  } else abundance != 0 | support != 0
  data.frame(
    taxon = as.character(feature_truth$taxon), direct_target = direct,
    truth_abundance_coefficient = abundance,
    truth_support_coefficient = support,
    truth_direction = sign(ifelse(abundance != 0, abundance, support)),
    check.names = FALSE
  )
}

.mb_evaluate_calls <- function(calls, feature_truth, applicability) {
  truth <- .mb_da_truth(feature_truth)
  table <- merge(calls, truth, by = "taxon", sort = FALSE)
  table <- table[match(calls$taxon, table$taxon), , drop = FALSE]
  significant <- table$significant
  target <- table$direct_target
  tp <- sum(significant & target)
  outside <- sum(significant & !target)
  fn <- sum(!significant & target)
  tn <- sum(!significant & !target)
  discoveries <- sum(significant)
  targets <- sum(target)
  nulls <- sum(!target)
  strict <- isTRUE(applicability$da_fdr)
  power_supported <- isTRUE(applicability$da_power) && targets > 0
  table$outcome <- ifelse(
    significant & target, "true_positive",
    ifelse(significant & !target, "outside_target_discovery",
           ifelse(!significant & target, "false_negative", "not_discovered_non_target"))
  )
  if (strict) {
    table$outcome[significant & !target] <- "false_positive"
    table$outcome[!significant & !target] <- "true_negative"
  }
  sign_accuracy <- NA_real_
  signed <- significant & target & table$truth_direction != 0
  if ("effect" %in% names(table) && any(signed)) {
    valid <- signed & !is.na(table$effect) & table$effect != 0
    if (any(valid)) sign_accuracy <- mean(sign(table$effect[valid]) == table$truth_direction[valid])
  }
  outside_fraction <- if (discoveries) outside / discoveries else 0
  metrics <- data.frame(
    da_truth_mode = as.character(applicability$da_truth_mode %||% "unsupported"),
    strict_non_target_null = strict, n_taxa = nrow(table),
    n_tested = sum(table$tested), coverage = mean(table$tested),
    n_direct_targets = targets, n_strict_null_taxa = if (strict) nulls else NA_integer_,
    n_significant = discoveries, direct_target_discoveries = tp,
    outside_target_discoveries = outside, missed_direct_targets = fn,
    outside_target_fraction_descriptive = outside_fraction,
    direct_target_power = if (power_supported) tp / targets else NA_real_,
    type_ii_error_rate = if (power_supported) fn / targets else NA_real_,
    direct_target_precision = if (discoveries) tp / discoveries else NA_real_,
    direct_target_sign_accuracy = sign_accuracy,
    true_positives = if (strict) tp else NA_integer_,
    false_positives = if (strict) outside else NA_integer_,
    true_negatives = if (strict) tn else NA_integer_,
    false_negatives = if (strict && targets) fn else NA_integer_,
    type_i_error_rate = if (strict && nulls) outside / nulls else NA_real_,
    specificity = if (strict && nulls) tn / nulls else NA_real_,
    fdp = if (strict) outside / max(discoveries, 1) else NA_real_,
    precision = if (strict && discoveries) tp / discoveries else NA_real_,
    fwer_indicator = if (strict) outside > 0 else NA,
    check.names = FALSE
  )
  list(metrics = metrics, table = table)
}

#' Build a complete DA decision table
#'
#' @param benchmark An mbBench benchmark object or bundle directory.
#' @param significant_taxa Character vector of taxa called significant.
#' @param tested_taxa Optional character vector of tested taxa. The default is all taxa.
#' @return A data frame with one row per benchmark taxon and logical `tested`
#'   and `significant` columns.
#' @export
mb_da_result <- function(benchmark, significant_taxa, tested_taxa = NULL) {
  bundle <- .mb_eval_bundle(benchmark)
  taxa <- as.character(bundle$feature_truth$taxon)
  significant_taxa <- as.character(significant_taxa)
  tested_taxa <- if (is.null(tested_taxa)) taxa else as.character(tested_taxa)
  unknown <- setdiff(union(significant_taxa, tested_taxa), taxa)
  if (length(unknown)) {
    .mb_abort("Unknown taxa: %s.", paste(utils::head(unknown, 5L), collapse = ", "))
  }
  if (length(setdiff(significant_taxa, tested_taxa))) {
    .mb_abort("Every significant taxon must also be tested.")
  }
  data.frame(taxon = taxa, tested = taxa %in% tested_taxa,
             significant = taxa %in% significant_taxa, check.names = FALSE)
}

#' Evaluate differential-abundance decisions
#'
#' `input` is a data frame with `taxon` and logical `significant` columns, a
#' named logical vector, a CSV path, or a named list of those forms. Optional
#' columns are `tested`, `effect`, `p_value`, and `q_value`. Every benchmark
#' taxon must be present; use [mb_da_result()] to construct a complete table.
#'
#' @param benchmark An mbBench benchmark object or bundle directory.
#' @param input One DA result or a named list of method results.
#' @param output_dir Optional directory for evaluation tables.
#' @param overwrite Replace a nonempty output directory.
#' @return An `mbbench_da_evaluation` list containing `metrics`,
#'   `taxon_results`, `config`, `output_dir`, and `files`.
#' @export
mb_evaluate_da <- function(benchmark, input, output_dir = NULL, overwrite = FALSE) {
  bundle <- .mb_eval_bundle(benchmark)
  taxa <- as.character(bundle$feature_truth$taxon)
  applicability <- bundle$config$truth_applicability %||% list()
  identity <- .mb_eval_identity(bundle)
  inputs <- .mb_da_inputs(input)
  metrics <- tables <- list()
  for (method in names(inputs)) {
    calls <- .mb_align_da(inputs[[method]], taxa, method)
    evaluated <- .mb_evaluate_calls(calls, bundle$feature_truth, applicability)
    metrics[[method]] <- cbind(
      data.frame(method = method, scenario = identity$scenario, seed = identity$seed),
      evaluated$metrics
    )
    table <- evaluated$table
    table$scenario <- identity$scenario
    table$seed <- identity$seed
    tables[[method]] <- table[c("method", "taxon", "scenario", "seed",
                                setdiff(names(table), c("method", "taxon", "scenario", "seed")))]
  }
  metrics <- .mb_bind_rows(metrics)
  tables <- .mb_bind_rows(tables)
  config <- list(schema_version = "mbbench_da_evaluation_v1",
                 benchmark = list(path = .mb_eval_bundle_path(bundle),
                                  scenario = identity$scenario, seed = identity$seed),
                 methods = names(inputs), truth_applicability = applicability,
                 full_taxon_truth_roster_used = TRUE)
  written <- .mb_eval_write(output_dir, overwrite, list(
    da_metrics.csv = metrics,
    da_taxon_results.csv.gz = tables,
    evaluation.yaml = config
  ))
  if (length(written$files)) {
    written$files <- written$files[c("da_metrics", "da_taxon_results", "evaluation")]
    names(written$files) <- c("metrics", "taxon_results", "config")
  }
  structure(list(metrics = metrics, taxon_results = tables,
                 applicability = applicability, config = config,
                 output_dir = written$output_dir, files = written$files,
                 summary = NULL),
            class = c("mbbench_da_evaluation", "list"))
}

.mb_metric_frames <- function(input, filename) {
  if (inherits(input, c("mbbench_da_evaluation", "mbbench_batch_evaluation"))) {
    return(input$metrics)
  }
  if (is.data.frame(input)) return(input)
  if (is.list(input) && !is.character(input)) {
    return(do.call(rbind, lapply(input, .mb_metric_frames, filename = filename)))
  }
  if (!is.character(input) || length(input) != 1L) {
    .mb_abort("input must be an evaluation, data frame, CSV, directory, or list.")
  }
  if (file.exists(input) && !dir.exists(input)) {
    if (identical(filename, "batch_metrics.csv") &&
        identical(basename(input), "batch_scores.csv")) {
      .mb_abort(paste0(
        "batch_scores.csv cannot be summarized without its batch_metrics.csv. ",
        "Pass the evaluation directory or batch_metrics.csv instead."
      ))
    }
    return(.mb_eval_frame(input))
  }
  files <- list.files(input, pattern = paste0("^", filename, "$"), recursive = TRUE,
                      full.names = TRUE)
  if (!length(files)) .mb_abort("No %s files were found below %s.", filename, input)
  do.call(rbind, lapply(files, .mb_eval_frame))
}

.mb_batch_score_frames <- function(input) {
  if (inherits(input, "mbbench_batch_evaluation")) return(input$scores)
  if (inherits(input, "mbbench_da_evaluation") || is.data.frame(input)) return(NULL)
  if (is.list(input) && !is.character(input)) {
    frames <- lapply(input, .mb_batch_score_frames)
    present <- !vapply(frames, is.null, logical(1))
    if (any(present) && !all(present)) {
      .mb_abort("Batch summary input mixes scored and unscored evaluations.")
    }
    return(if (any(present)) do.call(rbind, frames) else NULL)
  }
  if (!is.character(input) || length(input) != 1L) return(NULL)
  if (file.exists(input) && !dir.exists(input)) {
    candidate <- if (basename(input) == "batch_scores.csv") input else
      file.path(dirname(input), "batch_scores.csv")
    return(if (file.exists(candidate)) .mb_eval_frame(candidate) else NULL)
  }
  files <- list.files(input, pattern = "^batch_scores\\.csv$", recursive = TRUE,
                      full.names = TRUE)
  if (length(files)) do.call(rbind, lapply(files, .mb_eval_frame)) else NULL
}

.mb_bootstrap_interval <- function(values, repetitions) {
  values <- as.numeric(values[!is.na(values)])
  if (!length(values)) return(c(NA_real_, NA_real_))
  means <- replicate(repetitions, mean(sample(values, length(values), replace = TRUE)))
  as.numeric(stats::quantile(means, c(0.025, 0.975), names = FALSE))
}

.mb_bootstrap_seed <- function(seed, ...) {
  state <- as.double(seed) %% 2147483647
  bytes <- as.integer(charToRaw(enc2utf8(paste(..., sep = "\x1f"))))
  for (byte in bytes) {
    state <- (state * 1664525 + byte + 1013904223) %% 2147483647
  }
  as.integer(if (state <= 0) 1 else state)
}

.mb_wilson <- function(successes, total) {
  if (!total) return(c(NA_real_, NA_real_))
  z <- 1.959963984540054
  p <- successes / total
  denominator <- 1 + z ^ 2 / total
  center <- (p + z ^ 2 / (2 * total)) / denominator
  margin <- z * sqrt(p * (1 - p) / total + z ^ 2 / (4 * total ^ 2)) / denominator
  c(max(0, center - margin), min(1, center + margin))
}

.mb_summarize_da_frame <- function(metrics, bootstrap_repetitions, seed) {
  required <- c("method", "scenario", "seed", "da_truth_mode")
  missing <- setdiff(required, names(metrics))
  if (length(missing)) {
    .mb_abort("DA metrics are missing columns: %s.", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(metrics[c("method", "scenario", "seed")])) {
    .mb_abort("DA metrics repeat a method, scenario, and seed identity.")
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  groups <- split(seq_len(nrow(metrics)), interaction(metrics$method, metrics$scenario,
                                                   drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(index) {
    group <- metrics[index, , drop = FALSE]
    group <- group[order(as.character(group$seed), method = "radix"), , drop = FALSE]
    modes <- unique(as.character(stats::na.omit(group$da_truth_mode)))
    if (length(modes) != 1L) .mb_abort("DA truth mode varies within a summary group.")
    row <- list(method = group$method[[1L]], scenario = group$scenario[[1L]],
                n_repetitions = nrow(group), da_truth_mode = modes[[1L]])
    pairs <- c(fdp = "empirical_fdr", direct_target_power = "mean_power",
               type_i_error_rate = "mean_type_i_error_rate",
               type_ii_error_rate = "mean_type_ii_error_rate")
    for (source in names(pairs)) {
      target <- unname(pairs[[source]])
      values <- if (source %in% names(group)) group[[source]] else numeric()
      row[[target]] <- if (any(!is.na(values))) mean(values, na.rm = TRUE) else NA_real_
      set.seed(.mb_bootstrap_seed(seed, row$method, row$scenario, target))
      interval <- .mb_bootstrap_interval(values, bootstrap_repetitions)
      row[[paste0(target, "_ci_lower")]] <- interval[[1L]]
      row[[paste0(target, "_ci_upper")]] <- interval[[2L]]
    }
    valid <- if ("fwer_indicator" %in% names(group))
      group$fwer_indicator[!is.na(group$fwer_indicator)] else logical()
    row$fwer <- if (length(valid)) mean(valid) else NA_real_
    interval <- .mb_wilson(sum(valid), length(valid))
    row$fwer_ci_lower <- interval[[1L]]
    row$fwer_ci_upper <- interval[[2L]]
    as.data.frame(row, check.names = FALSE)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Summarize repeated DA evaluations
#'
#' Empirical FDR is the mean seed-level FDP. The function accepts an evaluation,
#' metrics data frame, CSV, directory containing `da_metrics.csv` files, or a list.
#'
#' @param input Evaluation results, a metrics table/path, directory, or list.
#' @param output_dir Optional summary directory.
#' @param overwrite Replace a nonempty output directory.
#' @param bootstrap_repetitions Number of deterministic bootstrap repetitions.
#' @param seed Bootstrap seed.
#' @return One summary row per method and scenario.
#' @export
mb_summarize_da <- function(input, output_dir = NULL, overwrite = FALSE,
                            bootstrap_repetitions = 2000L, seed = 42L) {
  bootstrap_repetitions <- .mb_scalar_number(
    bootstrap_repetitions, "bootstrap_repetitions", minimum = 1, integer = TRUE
  )
  seed <- .mb_scalar_number(seed, "seed", minimum = 0, integer = TRUE)
  metrics <- .mb_metric_frames(input, "da_metrics.csv")
  summary <- .mb_summarize_da_frame(metrics, bootstrap_repetitions, seed)
  .mb_eval_write(output_dir, overwrite, list(
    replicate_metrics.csv = metrics, summary.csv = summary,
    summary.yaml = list(schema_version = "mbbench_da_summary_v1",
                        bootstrap_repetitions = bootstrap_repetitions,
                        bootstrap_seed = seed)
  ))
  summary
}

#' Evaluate many DA benchmark repetitions from a manifest
#'
#' @param manifest Data frame or CSV with `benchmark`, `input`, and `method` columns.
#'   Relative paths in a CSV are resolved from the manifest directory.
#' @param output_dir Optional evaluation directory.
#' @param overwrite Replace a nonempty output directory.
#' @param bootstrap_repetitions Bootstrap repetitions for the included summary.
#' @param seed Bootstrap seed.
#' @return An `mbbench_da_evaluation` with repetition metrics and `summary`.
#' @export
mb_evaluate_da_many <- function(manifest, output_dir = NULL, overwrite = FALSE,
                                bootstrap_repetitions = 2000L, seed = 42L) {
  bootstrap_repetitions <- .mb_scalar_number(
    bootstrap_repetitions, "bootstrap_repetitions", minimum = 1, integer = TRUE
  )
  seed <- .mb_scalar_number(seed, "seed", minimum = 0, integer = TRUE)
  manifest_path <- if (is.character(manifest) && length(manifest) == 1L) {
    .mb_normalize_path(manifest, must_work = TRUE)
  } else NULL
  frame <- if (is.null(manifest_path)) manifest else .mb_eval_frame(manifest_path)
  if (!is.data.frame(frame) || !nrow(frame) ||
      !all(c("benchmark", "input", "method") %in% names(frame))) {
    .mb_abort("manifest must be a nonempty data frame with benchmark, input, and method columns.")
  }
  for (column in c("benchmark", "input", "method")) {
    values <- as.character(frame[[column]])
    if (anyNA(values) || any(!nzchar(trimws(values)))) {
      .mb_abort("manifest column '%s' cannot contain missing values.", column)
    }
    frame[[column]] <- values
  }
  if (!is.null(manifest_path)) {
    base <- dirname(manifest_path)
    for (column in c("benchmark", "input")) {
      frame[[column]] <- vapply(frame[[column]], function(path) {
        if (.Platform$file.sep == "\\" && grepl("^[A-Za-z]:[/\\\\]", path)) {
          .mb_normalize_path(path, must_work = TRUE)
        } else if (grepl("^[/\\\\]", path)) {
          .mb_normalize_path(path, must_work = TRUE)
        } else .mb_normalize_path(file.path(base, path), must_work = TRUE)
      }, character(1))
    }
  }
  evaluations <- lapply(seq_len(nrow(frame)), function(index) {
    named_input <- list(.mb_eval_frame(frame$input[[index]]))
    names(named_input) <- as.character(frame$method[[index]])
    mb_evaluate_da(frame$benchmark[[index]], named_input)
  })
  metrics <- do.call(rbind, lapply(evaluations, `[[`, "metrics"))
  taxon_results <- do.call(rbind, lapply(evaluations, `[[`, "taxon_results"))
  rownames(metrics) <- rownames(taxon_results) <- NULL
  summary <- .mb_summarize_da_frame(metrics, bootstrap_repetitions, seed)
  config <- list(
    schema_version = "mbbench_da_evaluation_collection_v1",
    n_runs = nrow(metrics), manifest = manifest_path,
    summary = list(bootstrap_repetitions = bootstrap_repetitions,
                   bootstrap_seed = seed)
  )
  written <- .mb_eval_write(output_dir, overwrite, list(
    da_metrics.csv = metrics, da_taxon_results.csv.gz = taxon_results,
    summary.csv = summary, evaluation.yaml = config
  ))
  if (length(written$files)) {
    written$files <- written$files[
      c("da_metrics", "da_taxon_results", "summary", "evaluation")
    ]
    names(written$files) <- c("metrics", "taxon_results", "summary", "config")
  }
  structure(list(metrics = metrics, taxon_results = taxon_results,
                 applicability = list(varies_by_benchmark = TRUE), config = config,
                 output_dir = written$output_dir, files = written$files,
                 summary = summary),
            class = c("mbbench_da_evaluation", "list"))
}

#' Summarize repeated batch evaluations
#'
#' @param input Evaluation results, a metrics table/path, directory, or list.
#' @param output_dir Optional summary directory.
#' @param overwrite Replace a nonempty output directory.
#' @return An `mbbench_batch_summary`. `replicate_metrics` and optional
#'   `replicate_scores` are wide tables with one row per method, scenario, and
#'   seed. `metric_summary` and optional `score_summary` are long tables with one
#'   row per method, scenario, and metric, including the repetition count, mean,
#'   and standard deviation. `replicate_scores` and `score_summary` are `NULL`
#'   when no `batch_scores.csv` files are present.
#' @export
mb_summarize_batch <- function(input, output_dir = NULL, overwrite = FALSE) {
  metrics <- .mb_metric_frames(input, "batch_metrics.csv")
  missing <- setdiff(c("method", "scenario", "seed"), names(metrics))
  if (length(missing)) {
    .mb_abort("Batch metrics are missing columns: %s.", paste(missing, collapse = ", "))
  }
  if (anyDuplicated(metrics[c("method", "scenario", "seed")])) {
    .mb_abort("Batch metrics repeat a method, scenario, and seed identity.")
  }
  summarize_numeric <- function(frame, identifiers) {
    numeric <- setdiff(names(frame)[vapply(frame, is.numeric, logical(1))], identifiers)
    groups <- split(seq_len(nrow(frame)), interaction(
      frame$method, frame$scenario, drop = TRUE, lex.order = TRUE
    ))
    rows <- list()
    for (index in groups) for (metric in numeric) {
      values <- frame[[metric]][index]
      values <- values[!is.na(values)]
      if (length(values)) rows[[length(rows) + 1L]] <- data.frame(
        method = frame$method[index[[1L]]],
        scenario = frame$scenario[index[[1L]]],
        metric = metric,
        n_repetitions = length(values),
        mean = mean(values),
        standard_deviation = if (length(values) > 1L) stats::sd(values) else NA_real_
      )
    }
    result <- if (length(rows)) do.call(rbind, rows) else data.frame()
    rownames(result) <- NULL
    result
  }
  metric_identifiers <- c(
    "method", "scenario", "seed", "batch_recovery_applicable",
    "biological_preservation_applicable", "relative_pseudocount",
    "structure_prevalence_min", "structure_maximum_relative_abundance_min",
    "structure_n_active_taxa", "structure_metrics_supported",
    "structure_relative_pseudocount", "structure_pca_components"
  )
  metric_summary <- summarize_numeric(metrics, metric_identifiers)
  scores <- .mb_batch_score_frames(input)
  score_summary <- NULL
  if (!is.null(scores)) {
    score_missing <- setdiff(c("method", "scenario", "seed"), names(scores))
    if (length(score_missing)) {
      .mb_abort("Batch scores are missing columns: %s.",
                paste(score_missing, collapse = ", "))
    }
    if (anyDuplicated(scores[c("method", "scenario", "seed")])) {
      .mb_abort("Batch scores repeat a method, scenario, and seed identity.")
    }
    identities <- c("method", "scenario", "seed")
    key <- function(frame) do.call(paste, c(lapply(frame[identities], as.character), sep = "\r"))
    if (!setequal(key(metrics), key(scores))) {
      .mb_abort("Batch summary input mixes scored and unscored evaluation identities.")
    }
    score_summary <- summarize_numeric(scores, c(
      "method", "scenario", "seed", "batch_recovery_applicable",
      "biological_preservation_applicable", "score_definition",
      "classifier_seed", "score_supported"
    ))
  }
  output_files <- list(
    replicate_metrics.csv = metrics, summary.csv = metric_summary,
    summary.yaml = list(schema_version = "mbbench_batch_summary_v1")
  )
  if (!is.null(scores)) {
    output_files$replicate_scores.csv <- scores
    output_files$score_summary.csv <- score_summary
  }
  written <- .mb_eval_write(output_dir, overwrite, output_files)
  files <- list()
  if (!is.null(written$output_dir)) {
    files <- list(
      replicate_metrics = file.path(written$output_dir, "replicate_metrics.csv"),
      metric_summary = file.path(written$output_dir, "summary.csv"),
      config = file.path(written$output_dir, "summary.yaml")
    )
    if (!is.null(scores)) {
      files$replicate_scores <- file.path(written$output_dir, "replicate_scores.csv")
      files$score_summary <- file.path(written$output_dir, "score_summary.csv")
    }
  }
  structure(list(
    metric_summary = metric_summary,
    score_summary = score_summary,
    replicate_metrics = metrics,
    replicate_scores = scores,
    output_dir = written$output_dir,
    files = files
  ), class = c("mbbench_batch_summary", "list"))
}
