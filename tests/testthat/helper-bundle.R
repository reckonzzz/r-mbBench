write_gz_csv <- function(data, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(data, con, row.names = FALSE, quote = TRUE)
}

make_test_bundle <- function(root = tempfile("mbbench-bundle-")) {
  dir.create(root, recursive = TRUE)
  counts <- data.frame(sample_id = c("sample_1", "sample_2"), taxon_a = c(2L, 3L), taxon_b = c(8L, 7L), check.names = FALSE)
  relative <- data.frame(sample_id = counts$sample_id, taxon_a = c(0.2, 0.3), taxon_b = c(0.8, 0.7), check.names = FALSE)
  latent <- data.frame(sample_id = counts$sample_id, taxon_a = c(0.25, 0.35), taxon_b = c(0.75, 0.65), check.names = FALSE)
  reference <- data.frame(sample_id = counts$sample_id, taxon_a = c(0.22, 0.32), taxon_b = c(0.78, 0.68), check.names = FALSE)
  metadata <- data.frame(sample_id = counts$sample_id, group = c("control", "case"), batch = c("batch_1", "batch_2"))
  feature_truth <- data.frame(
    taxon = c("taxon_a", "taxon_b"),
    is_direct_biological_target = c(TRUE, FALSE),
    biological_abundance_coefficient = c(0.7, 0),
    biological_support_coefficient = c(0, 0),
    is_batch_target = c(FALSE, TRUE)
  )
  sample_truth <- data.frame(sample_id = counts$sample_id, biological_effect = c(0, 1), depth = c(10L, 10L))
  tables <- list(
    counts.csv.gz = counts,
    relative_abundance.csv.gz = relative,
    latent_composition.csv.gz = latent,
    metadata.csv.gz = metadata,
    feature_truth.csv.gz = feature_truth,
    sample_truth.csv.gz = sample_truth,
    batch_neutral_biological_reference.csv.gz = reference
  )
  for (name in names(tables)) write_gz_csv(tables[[name]], file.path(root, name))
  yaml::write_yaml(list(
    schema_version = "mbbench_benchmark_v1", preset = "test",
    runtime = list(seed = 42L), design = list(n_samples = 2L),
    truth_applicability = list(
      da_power = TRUE, da_fdr = TRUE,
      da_truth_mode = "strict_relative_abundance_targets",
      batch_recovery = TRUE, biological_preservation = TRUE
    )
  ),
                   file.path(root, "benchmark.yaml"))
  root
}

zip_test_bundle <- function(bundle) {
  archive <- tempfile(fileext = ".zip")
  old <- setwd(dirname(bundle))
  on.exit(setwd(old), add = TRUE)
  status <- utils::zip(archive, files = basename(bundle), flags = "-rq")
  if (!identical(status, 0L) || !file.exists(archive)) return(NULL)
  archive
}
