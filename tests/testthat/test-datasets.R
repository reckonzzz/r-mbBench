test_that("an extracted benchmark bundle returns ordinary R tables", {
  bundle <- make_test_bundle()
  dat <- mb_dataset(path = bundle)

  expect_s3_class(dat, "mbbench_dataset")
  expect_s3_class(dat, "list")
  expect_equal(nrow(dat$counts), 2L)
  expect_identical(names(dat$counts), c("sample_id", "taxon_a", "taxon_b"))
  expect_equal(dat$relative_abundance$taxon_a, c(0.2, 0.3))
  expect_null(dat$batch_neutral_counts)
  expect_identical(dat$source$type, "local_directory")
  expect_output(print(dat), "2 samples x 2 taxa")
})

test_that("dataset source is unambiguous", {
  bundle <- make_test_bundle()
  expect_error(mb_dataset(), "exactly one")
  expect_error(mb_dataset(path = bundle, archive = "anything.zip"), "exactly one")
})

test_that("bundle validation reports structural and numeric errors", {
  bundle <- make_test_bundle()
  unlink(file.path(bundle, "feature_truth.csv.gz"))
  expect_error(mb_dataset(path = bundle), "missing: feature_truth.csv.gz")

  bundle <- make_test_bundle()
  bad <- data.frame(sample_id = c("sample_1", "sample_2"), taxon_a = c(0.1, 0.3), taxon_b = c(0.9, 0.7))
  write_gz_csv(bad, file.path(bundle, "relative_abundance.csv.gz"))
  expect_error(mb_dataset(path = bundle), "counts divided by row depth")
})

test_that("a manually downloaded ZIP can be imported without Python", {
  bundle <- make_test_bundle()
  archive <- zip_test_bundle(bundle)
  skip_if(is.null(archive), "No ZIP command is available")
  cache <- tempfile("mbbench-cache-")
  old <- options(mbBench.cache_dir = cache)
  on.exit(options(old), add = TRUE)

  dat <- mb_dataset(archive = archive)
  expect_s3_class(dat, "mbbench_dataset")
  expect_identical(dat$source$type, "manual_archive")
  expect_true(dir.exists(dat$path))
})

test_that("catalog schema rejects ambiguous releases", {
  valid <- list(
    schema_version = "mbbench_dataset_catalog_v1",
    latest = "1",
    releases = list(list(
      version = "1",
      datasets = list(list(
        preset = "da-null-balanced", seed = 42L,
        archive = "dataset.zip", url = "https://example.org/dataset.zip"
      ))
    ))
  )
  expect_invisible(mbBench:::.mb_validate_catalog(valid))
  duplicate <- valid
  duplicate$releases[[1]]$datasets[[2]] <- duplicate$releases[[1]]$datasets[[1]]
  expect_error(mbBench:::.mb_validate_catalog(duplicate), "repeats preset")
})
