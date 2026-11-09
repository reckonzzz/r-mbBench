make_evaluation_bundle <- function() {
  taxa <- c("taxon_a", "taxon_b", "taxon_c", "rare_taxon")
  design <- expand.grid(replicate = seq_len(4L), batch = c("batch_1", "batch_2"),
                        group = c("control", "case"), stringsAsFactors = FALSE)
  counts_matrix <- t(vapply(seq_len(nrow(design)), function(index) {
    group <- as.integer(design$group[[index]] == "case")
    batch <- as.integer(design$batch[[index]] == "batch_2")
    replicate <- design$replicate[[index]]
    c(40 + 7 * group + replicate, 25 + 9 * batch + replicate,
      30 + 3 * group + 2 * batch, as.integer(index == 1L))
  }, numeric(4)))
  sample_ids <- sprintf("sample_%02d", seq_len(nrow(design)))
  counts <- data.frame(sample_id = sample_ids, counts_matrix, check.names = FALSE)
  names(counts)[-1L] <- taxa
  relative_matrix <- counts_matrix / rowSums(counts_matrix)
  relative <- data.frame(sample_id = sample_ids, relative_matrix, check.names = FALSE)
  names(relative)[-1L] <- taxa
  reference_matrix <- relative_matrix
  reference_matrix[, 1L] <- reference_matrix[, 1L] + 0.03
  reference_matrix[, 2L] <- reference_matrix[, 2L] - 0.02
  reference_matrix[, 3L] <- reference_matrix[, 3L] - 0.01
  reference <- data.frame(sample_id = sample_ids, reference_matrix, check.names = FALSE)
  names(reference)[-1L] <- taxa
  metadata <- data.frame(sample_id = sample_ids, group = design$group,
                         batch = design$batch, stringsAsFactors = FALSE)
  feature_truth <- data.frame(
    taxon = taxa, is_direct_biological_target = c(TRUE, FALSE, TRUE, FALSE),
    biological_abundance_coefficient = c(0.7, 0, -0.7, 0),
    biological_support_coefficient = 0
  )
  list(
    counts = counts, relative_abundance = relative,
    latent_composition = relative, metadata = metadata,
    feature_truth = feature_truth,
    sample_truth = data.frame(sample_id = sample_ids),
    batch_neutral_reference = reference,
    config = list(
      preset = "strict", runtime = list(seed = 42L),
      truth_applicability = list(
        da_power = TRUE, da_fdr = TRUE,
        da_truth_mode = "strict_relative_abundance_targets",
        batch_recovery = TRUE, biological_preservation = TRUE
      )
    )
  )
}

test_that("batch evaluation preserves paired and structure metric boundaries", {
  bundle <- make_evaluation_bundle()
  identity <- mb_evaluate_batch(bundle, list(none = bundle$relative_abundance))
  expect_s3_class(identity, "mbbench_batch_evaluation")
  expect_null(identity$scores)
  expect_equal(identity$metrics$pbrs_raw, 0, tolerance = 1e-12)
  expect_equal(identity$metrics$fraction_samples_unchanged, 1)
  expect_equal(identity$metrics$structure_n_active_taxa, 3)
  expect_equal(identity$metrics$batch_asw_change, 0, tolerance = 1e-12)
  expect_true(identity$metrics$batch_recovery_applicable)

  reference <- mb_evaluate_batch(bundle, list(reference = bundle$batch_neutral_reference))
  expect_equal(reference$metrics$pbrs_raw, 1, tolerance = 1e-12)
  expect_equal(reference$metrics$corrected_rms_aitchison, 0, tolerance = 1e-12)

  rare_changed <- bundle$relative_abundance
  rare_changed$rare_taxon <- rare_changed$rare_taxon + 0.02
  changed <- mb_evaluate_batch(bundle, rare_changed)
  expect_equal(changed$metrics$structure_n_active_taxa, 3)
  expect_gt(changed$metrics$corrected_rms_aitchison, 0)

  excluded_only <- bundle$relative_abundance
  excluded_only[c("taxon_a", "taxon_b", "taxon_c")] <- 0
  excluded_only$rare_taxon <- 1
  mixed <- mb_evaluate_batch(
    bundle,
    list(identity = bundle$relative_abundance, excluded_only = excluded_only)
  )
  unsupported <- mixed$metrics[mixed$metrics$method == "excluded_only", ]
  expect_false(unsupported$structure_metrics_supported)
  expect_match(unsupported$structure_metrics_reason, "zero-mass sample")
})

test_that("batch scores are optional and require lisi only when selected", {
  bundle <- make_evaluation_bundle()
  expect_silent(mb_evaluate_batch(bundle, bundle$relative_abundance))
  if (!requireNamespace("lisi", quietly = TRUE)) {
    expect_error(
      mb_evaluate_batch(
        bundle, bundle$relative_abundance, calculate_scores = TRUE
      ),
      "requires the R package 'lisi'"
    )
  }
  expect_error(
    mb_evaluate_batch(bundle, bundle$relative_abundance, calculate_scores = NA),
    "must be TRUE or FALSE"
  )
})

test_that("R publication score transformations match declared formulas", {
  expect_equal(mbBench:::.mb_removal_fraction(0.2, 0.1), 0.5)
  expect_equal(mbBench:::.mb_mixing_improvement(0.4, 0.7), 0.5)
  expect_equal(mbBench:::.mb_bounded_similarity(0.4, 0.3), 0.75)
  expect_true(is.na(mbBench:::.mb_bounded_similarity(NULL, 0.3)))
})

test_that("DA evaluation uses the complete truth roster", {
  bundle <- make_evaluation_bundle()
  calls <- mb_da_result(bundle, c("taxon_a", "taxon_c"))
  result <- mb_evaluate_da(bundle, list(perfect = calls))
  expect_s3_class(result, "mbbench_da_evaluation")
  expect_equal(result$metrics$direct_target_power, 1)
  expect_equal(result$metrics$type_i_error_rate, 0)
  expect_equal(result$metrics$fdp, 0)
  expect_error(mb_evaluate_da(bundle, calls[-1L, ]), "do not match")
})

test_that("many DA repetitions can be summarized without a long list", {
  metrics <- data.frame(
    method = "method_a", scenario = "strict", seed = seq_len(100L),
    da_truth_mode = "strict_relative_abundance_targets",
    fdp = c(rep(0.2, 25L), rep(0, 75L)), direct_target_power = 0.5,
    type_i_error_rate = 0.01, type_ii_error_rate = 0.5,
    fwer_indicator = c(rep(TRUE, 10L), rep(FALSE, 90L))
  )
  summary <- mb_summarize_da(metrics, bootstrap_repetitions = 100L)
  expect_equal(summary$empirical_fdr, 0.05)
  expect_equal(summary$mean_power, 0.5)
  expect_equal(summary$fwer, 0.1)
})

test_that("DA bootstrap intervals are stable across methods and row order", {
  rows <- function(method, values) data.frame(
    method = method, scenario = "strict", seed = seq_along(values) - 1L,
    da_truth_mode = "strict_relative_abundance_targets",
    fdp = values, direct_target_power = 1 - values,
    type_i_error_rate = values / 2, type_ii_error_rate = values / 3,
    fwer_indicator = values > 0.25
  )
  target <- rows("b", c(0, 0.1, 0.4, 0.8, 0.2))
  other <- rows("a", c(0.9, 0.8, 0.7, 0.6, 0.5))
  alone <- mb_summarize_da(target, bootstrap_repetitions = 200L)
  together <- mb_summarize_da(rbind(other, target), bootstrap_repetitions = 200L)
  together <- together[together$method == "b", , drop = FALSE]
  reversed <- mb_summarize_da(target[nrow(target):1L, ], bootstrap_repetitions = 200L)
  intervals <- grepl("_ci_(lower|upper)$", names(alone))
  expect_equal(unlist(alone[intervals], use.names = FALSE),
               unlist(together[intervals], use.names = FALSE))
  expect_equal(unlist(alone[intervals], use.names = FALSE),
               unlist(reversed[intervals], use.names = FALSE))
})

test_that("batch summaries keep raw metrics and scores separate", {
  root <- tempfile("mbbench_batch_summary_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  for (seed in c(42L, 43L)) {
    run <- file.path(root, as.character(seed))
    dir.create(run)
    utils::write.csv(data.frame(
      method = "method_a", scenario = "scenario_a", seed = seed,
      pbrs = 0.4 + 0.1 * (seed - 42L)
    ), file.path(run, "batch_metrics.csv"), row.names = FALSE)
    utils::write.csv(data.frame(
      method = "method_a", scenario = "scenario_a", seed = seed,
      overall_score = 0.6 + 0.1 * (seed - 42L)
    ), file.path(run, "batch_scores.csv"), row.names = FALSE)
  }
  summary <- mb_summarize_batch(root)
  expect_s3_class(summary, "mbbench_batch_summary")
  expect_equal(summary$metric_summary$mean[summary$metric_summary$metric == "pbrs"], 0.45)
  expect_equal(summary$score_summary$mean[summary$score_summary$metric == "overall_score"], 0.65)
})

test_that("batch summaries reject partly scored repetitions", {
  root <- tempfile("mbbench_partial_scores_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  for (seed in c(42L, 43L)) {
    run <- file.path(root, as.character(seed))
    dir.create(run)
    utils::write.csv(data.frame(
      method = "method_a", scenario = "scenario_a", seed = seed
    ), file.path(run, "batch_metrics.csv"), row.names = FALSE)
    if (seed == 42L) utils::write.csv(data.frame(
      method = "method_a", scenario = "scenario_a", seed = seed,
      overall_score = 0.6
    ), file.path(run, "batch_scores.csv"), row.names = FALSE)
  }
  expect_error(mb_summarize_batch(root), "mixes scored and unscored")
})

test_that("batch score files cannot be used as metric inputs", {
  root <- tempfile("batch_score_input_")
  dir.create(root)
  path <- file.path(root, "batch_scores.csv")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  utils::write.csv(data.frame(
    method = "method_a", scenario = "scenario_a", seed = 42L,
    overall_score = 0.7
  ), path, row.names = FALSE)
  expect_error(mb_summarize_batch(path), "cannot be summarized without")
})

test_that("evaluation output replacement is staged and overwrite is strict", {
  output <- tempfile("mbbench_eval_output_")
  dir.create(output)
  sentinel <- file.path(output, "sentinel.txt")
  writeLines("keep", sentinel)
  on.exit(unlink(output, recursive = TRUE, force = TRUE), add = TRUE)
  expect_error(
    mbBench:::.mb_eval_write(output, "False", list(summary.csv = data.frame(x = 1))),
    "overwrite must be TRUE or FALSE"
  )
  expect_error(
    mbBench:::.mb_eval_write(
      output, TRUE,
      list(first.csv = data.frame(x = 1), bad.csv = new.env())
    )
  )
  expect_equal(readLines(sentinel), "keep")

  file_output <- tempfile("mbbench_eval_file_")
  writeLines("keep", file_output)
  on.exit(unlink(file_output, force = TRUE), add = TRUE)
  expect_error(
    mbBench:::.mb_eval_write(
      file_output, FALSE, list(summary.csv = data.frame(x = 1))
    ),
    "existing file"
  )
  expect_equal(readLines(file_output), "keep")
})

test_that("bundle validation rejects duplicate sample IDs", {
  bundle <- make_evaluation_bundle()
  for (name in c("counts", "relative_abundance", "latent_composition", "metadata",
                 "sample_truth", "batch_neutral_reference")) {
    bundle[[name]]$sample_id[[2L]] <- bundle[[name]]$sample_id[[1L]]
  }
  expect_error(
    mbBench:::.mb_validate_bundle(bundle, bundle$config),
    "sample IDs must be nonempty and unique"
  )
})
