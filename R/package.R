#' mbBench: conditional microbiome simulation from R
#'
#' Pre-generated benchmark datasets can be downloaded, loaded, and evaluated
#' entirely in R. New data generation uses the optional mbBench Python backend
#' through `reticulate`.
#'
#' Start with [mb_dataset()] for pre-generated data, [mb_load_model()] for a
#' model bundle, [mb_generate()] for ordinary conditional generation, or
#' [mb_benchmark()] for truth-informed benchmark generation,
#' [mb_evaluate_batch()] for batch correction, or [mb_evaluate_da()] for DA.
#'
#' @keywords internal
"_PACKAGE"

.mb_backend_requirement <- "mbBench>=0.1.0,<0.2.0"
.mb_catalog_schema <- "mbbench_dataset_catalog_v1"
.mb_catalog_filename <- "mbbench_dataset_catalog.json"
.mb_official_dataset_record_url <- "https://zenodo.org/api/records/24076644"
