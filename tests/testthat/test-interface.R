test_that("backend status does not initialize Python by default", {
  cache <- tempfile("mbbench-cache-")
  old <- options(mbBench.cache_dir = cache)
  on.exit(options(old), add = TRUE)

  status <- mb_backend_status()
  expect_type(status, "list")
  expect_false(status$configured)
  expect_false(status$checked)
  expect_true(is.na(status$available))
  expect_named(status, c("configured", "initialized", "checked", "available", "requirement", "python", "version"))
})

test_that("model and scenario inputs fail early with clear messages", {
  expect_error(mb_generate(), "exactly one of model or bank")
  expect_error(mb_generate(list()), "mb_load_model")
  fake_model <- structure(list(simulator = "placeholder"), class = c("mbbench_model", "list"))
  fake_bank <- structure(list(bank = "placeholder"), class = c("mbbench_bank", "list"))
  expect_error(mb_generate(fake_model, bank = fake_bank), "exactly one of model or bank")
  expect_error(mb_generate(bank = list()), "mb_load_bank")
  expect_error(mb_benchmark(fake_model), "exactly one of preset or config")
  expect_error(
    mb_benchmark(fake_model, preset = "x", config = "x.yaml"),
    "exactly one of preset or config"
  )
  expect_error(
    mb_generate(fake_model, load = FALSE),
    "load = FALSE requires output_dir"
  )
  expect_error(
    mb_generate_from_condition_table(fake_model, data.frame()),
    "nonempty data.frame"
  )
  expect_error(
    mb_benchmark(fake_model, preset = "x", load = FALSE),
    "load = FALSE requires output_dir"
  )
})
