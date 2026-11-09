.mb_python_value <- function(value) {
  converted <- reticulate::py_to_r(value)
  converted
}

.mb_python_path <- function(value) {
  converted <- reticulate::py_to_r(value)
  if (is.null(converted)) return(NULL)
  if (is.character(converted)) return(as.character(converted))
  as.character(reticulate::py_to_r(value$`__str__`()))
}

.mb_python_files <- function(value) {
  files <- reticulate::py_to_r(value)
  if (is.null(files) || !length(files)) return(stats::setNames(character(), character()))
  vapply(files, .mb_python_path, character(1))
}

.mb_generation_output <- function(result) {
  abundance <- .mb_python_value(result$abundance)
  counts <- .mb_python_value(result$counts)
  metadata <- .mb_python_value(result$metadata)
  structure(list(
    abundance = abundance,
    counts = counts,
    metadata = metadata,
    config = .mb_python_value(result$config),
    files = .mb_python_files(result$files),
    output_dir = .mb_python_path(result$output_dir),
    loaded = !is.null(abundance),
    estimated_memory_mb = as.numeric(.mb_python_value(result$estimated_memory_mb))
  ), class = c("mbbench_generation", "list"))
}

.mb_generation_settings <- function(output_mode, guidance_scale, support_temperature,
                                    seed, library_size, use_body_system, batch_size,
                                    output_dir, overwrite, load, max_in_memory_mb) {
  output_mode <- match.arg(output_mode, c("relative", "observed"))
  guidance_scale <- .mb_scalar_number(guidance_scale, "guidance_scale", minimum = 0)
  support_temperature <- .mb_scalar_number(support_temperature, "support_temperature", minimum = 0)
  if (support_temperature <= 0) .mb_abort("support_temperature must be positive.")
  seed <- .mb_scalar_number(seed, "seed", minimum = 0, integer = TRUE)
  if (!is.null(library_size)) {
    library_size <- .mb_scalar_number(library_size, "library_size", minimum = 1, integer = TRUE)
  }
  if (!is.null(use_body_system) && (!is.logical(use_body_system) || length(use_body_system) != 1L || is.na(use_body_system))) {
    .mb_abort("use_body_system must be TRUE, FALSE, or NULL.")
  }
  batch_size <- .mb_scalar_number(batch_size, "batch_size", minimum = 1, integer = TRUE)
  max_in_memory_mb <- .mb_scalar_number(max_in_memory_mb, "max_in_memory_mb", minimum = 0)
  if (max_in_memory_mb <= 0) .mb_abort("max_in_memory_mb must be positive.")
  if (!is.null(output_dir)) output_dir <- .mb_normalize_path(output_dir, must_work = FALSE)
  if (!is.logical(load) || length(load) != 1L || is.na(load)) {
    .mb_abort("load must be TRUE or FALSE.")
  }
  overwrite <- .mb_scalar_logical(overwrite, "overwrite")
  if (!load && is.null(output_dir)) {
    .mb_abort("load = FALSE requires output_dir for streamed generation.")
  }
  list(
    output_mode = output_mode,
    guidance_scale = guidance_scale,
    support_temperature = support_temperature,
    seed = seed,
    library_size = library_size,
    use_body_system = use_body_system,
    batch_size = batch_size,
    output_dir = output_dir,
    overwrite = overwrite,
    load = load,
    max_in_memory_mb = max_in_memory_mb
  )
}

#' Generate conditional microbiome samples
#'
#' Uses either a model returned by [mb_load_model()] or a bank returned by
#' [mb_load_bank()]. Small results are returned as R
#' data frames. If the Python backend estimates that the request exceeds
#' `max_in_memory_mb`, provide `output_dir` for streamed writing or raise the
#' limit when enough memory is available.
#'
#' @param model An `mbbench_model` object, or `NULL` when `bank` is used.
#' @param bank An `mbbench_bank` object, or `NULL` when `model` is used.
#' @param phenotype Phenotype condition, used only when the selected model or
#'   bank contains the standard `Phenotype_Tokens` field.
#' @param body_site Body-site condition, used only when the selected model or
#'   bank contains the standard `Body_Site_Coarse` field.
#' @param n_samples Number of samples.
#' @param conditions Named list of any additional model condition fields.
#' @param output_mode `"relative"` for relative abundance or `"observed"` for
#'   observed counts plus relative abundance.
#' @param guidance_scale Classifier-free guidance scale.
#' @param support_temperature Support sampling temperature.
#' @param seed Random seed.
#' @param library_size Optional fixed library size for observed output.
#' @param use_body_system Whether to use the model's derived Body-System
#'   condition. `NULL` lets the backend decide from the model.
#' @param batch_size Generation chunk size.
#' @param output_dir Optional directory for generated CSV files and the resolved
#'   `generation.yaml`.
#' @param overwrite Replace a nonempty output directory.
#' @param load Whether to return generated tables in R when they fit within
#'   `max_in_memory_mb`. If `FALSE`, `output_dir` is required and results are
#'   only written to files. Large requests are streamed even when `TRUE`.
#' @param max_in_memory_mb Maximum estimated memory for returning full tables.
#' @return An `mbbench_generation` list containing `abundance`, optional
#'   `counts`, `metadata`, `config`, `files`, `output_dir`, and
#'   `estimated_memory_mb`. `loaded` reports whether full tables were returned;
#'   streamed runs leave those table fields as `NULL`.
#' @export
#' @examples
#' \dontrun{
#' model <- mb_load_model(path = "my_model_bundle")
#' generated <- mb_generate(model, n_samples = 20)
#' generated$abundance
#'
#' mb_generate(model, n_samples = 100000, output_dir = "large_run")
#' }
mb_generate <- function(model = NULL, bank = NULL,
                        phenotype = "Healthy", body_site = "stool_feces",
                        n_samples = 100L, conditions = list(),
                        output_mode = c("relative", "observed"),
                        guidance_scale = 1.25, support_temperature = 0.60,
                        seed = 42L, library_size = NULL, use_body_system = NULL,
                        batch_size = 64L, output_dir = NULL, overwrite = FALSE,
                        load = TRUE,
                        max_in_memory_mb = 512) {
  source <- .mb_select_source(model, bank)
  phenotype <- .mb_scalar_string(phenotype, "phenotype")
  body_site <- .mb_scalar_string(body_site, "body_site")
  n_samples <- .mb_scalar_number(n_samples, "n_samples", minimum = 1, integer = TRUE)
  if (length(conditions)) conditions <- .mb_named_list(conditions, "conditions")
  settings <- .mb_generation_settings(
    output_mode, guidance_scale, support_temperature, seed, library_size,
    use_body_system, batch_size, output_dir, overwrite, load, max_in_memory_mb
  )
  backend <- .mb_ensure_backend()
  arguments <- list(
    phenotype = phenotype, body_site = body_site, n_samples = n_samples,
    conditions = conditions, output_mode = settings$output_mode,
    guidance_scale = settings$guidance_scale,
    support_temperature = settings$support_temperature, seed = settings$seed,
    library_size = settings$library_size, use_body_system = settings$use_body_system,
    batch_size = settings$batch_size, output_dir = settings$output_dir,
    overwrite = settings$overwrite, load = settings$load,
    max_in_memory_mb = settings$max_in_memory_mb
  )
  if (source$type == "model") arguments$model <- source$value$simulator
  else arguments$bank <- source$value$bank
  result <- do.call(backend$generate, arguments)
  .mb_generation_output(result)
}

#' Generate samples from a condition table
#'
#' Each row of `condition_table` describes one condition and may contain
#' `n_samples` plus one column per model condition. The `phenotype` and
#' `body_site` aliases are accepted only when the corresponding model fields
#' exist. Results use the same structure and memory rules as [mb_generate()].
#'
#' @param model An `mbbench_model` object, or `NULL` when `bank` is used.
#' @param condition_table A nonempty data frame with one condition combination
#'   per row. Use one column per model condition and optional `n_samples`.
#'   Multilabel cells may use `"A|B"` or `"A;B"`; list columns may contain
#'   character vectors.
#' @inheritParams mb_generate
#' @return An `mbbench_generation` list.
#' @export
#' @examples
#' \dontrun{
#' model <- mb_load_model(path = "my_model_bundle")
#' condition_table <- data.frame(
#'   phenotype = c("Healthy", "Disease"),
#'   body_site = "stool_feces",
#'   n_samples = c(20L, 20L)
#' )
#' generated <- mb_generate_from_condition_table(model, condition_table)
#' }
mb_generate_from_condition_table <- function(model = NULL, condition_table,
                                             bank = NULL,
                                             output_mode = c("relative", "observed"),
                                             guidance_scale = 1.25,
                                             support_temperature = 0.60,
                                             seed = 42L, library_size = NULL,
                                             use_body_system = NULL,
                                             batch_size = 64L, output_dir = NULL,
                                             overwrite = FALSE, load = TRUE,
                                             max_in_memory_mb = 512) {
  source <- .mb_select_source(model, bank)
  if (!is.data.frame(condition_table) || !nrow(condition_table)) {
    .mb_abort("condition_table must be a nonempty data.frame.")
  }
  duplicates <- unique(names(condition_table)[duplicated(names(condition_table))])
  if (length(duplicates)) {
    .mb_abort(
      "condition_table contains duplicate columns: %s.",
      paste(duplicates, collapse = ", ")
    )
  }
  settings <- .mb_generation_settings(
    output_mode, guidance_scale, support_temperature, seed, library_size,
    use_body_system, batch_size, output_dir, overwrite, load, max_in_memory_mb
  )
  backend <- .mb_ensure_backend()
  arguments <- list(
    condition_table = condition_table, output_mode = settings$output_mode,
    guidance_scale = settings$guidance_scale,
    support_temperature = settings$support_temperature, seed = settings$seed,
    library_size = settings$library_size, use_body_system = settings$use_body_system,
    batch_size = settings$batch_size, output_dir = settings$output_dir,
    overwrite = settings$overwrite, load = settings$load,
    max_in_memory_mb = settings$max_in_memory_mb
  )
  if (source$type == "model") arguments$model <- source$value$simulator
  else arguments$bank <- source$value$bank
  result <- do.call(backend$generate_from_condition_table, arguments)
  .mb_generation_output(result)
}

.mb_benchmark_output <- function(result) {
  counts <- .mb_python_value(result$counts)
  structure(list(
    counts = counts,
    relative_abundance = .mb_python_value(result$relative_abundance),
    latent_composition = .mb_python_value(result$latent_composition),
    metadata = .mb_python_value(result$metadata),
    feature_truth = .mb_python_value(result$feature_truth),
    sample_truth = .mb_python_value(result$sample_truth),
    batch_neutral_reference = .mb_python_value(result$batch_neutral_reference),
    batch_neutral_counts = .mb_python_value(result$batch_neutral_counts),
    config = .mb_python_value(result$config),
    files = .mb_python_files(result$files),
    output_dir = .mb_python_path(result$output_dir),
    loaded = !is.null(counts),
    estimated_memory_mb = as.numeric(.mb_python_value(result$estimated_memory_mb))
  ), class = c("mbbench_benchmark", "list"))
}

#' List benchmark presets available in the Python backend
#'
#' @return A data frame with `preset` and `description` columns.
#' @export
#' @examples
#' \dontrun{mb_list_presets()}
mb_list_presets <- function() {
  backend <- .mb_ensure_backend()
  rows <- reticulate::py_to_r(backend$list_presets())
  if (!length(rows)) {
    return(data.frame(preset = character(), description = character()))
  }
  data.frame(
    preset = vapply(rows, function(row) as.character(row[[1L]]), character(1)),
    description = vapply(rows, function(row) as.character(row[[2L]]), character(1)),
    stringsAsFactors = FALSE
  )
}

#' Show one benchmark preset
#'
#' @param preset Built-in benchmark preset name.
#' @return A named list containing the complete preset configuration.
#' @export
#' @examples
#' \dontrun{mb_show_preset("da-strict-common-moderate")}
mb_show_preset <- function(preset) {
  preset <- .mb_scalar_string(preset, "preset")
  backend <- .mb_ensure_backend()
  reticulate::py_to_r(backend$show_preset(preset))
}

#' List benchmark parameters
#'
#' Returns the same parameter names, defaults, accepted values, descriptions,
#' and YAML paths shown by `mbbench benchmark --list-parameters`.
#'
#' @param preset Optional preset whose defaults should be shown.
#' @return A data frame with one row per configurable parameter.
#' @export
#' @examples
#' \dontrun{mb_list_parameters("da-strict-common-moderate")}
mb_list_parameters <- function(preset = NULL) {
  if (!is.null(preset)) preset <- .mb_scalar_string(preset, "preset")
  backend <- .mb_ensure_backend()
  rows <- reticulate::py_to_r(backend$list_parameters(preset))
  if (!length(rows)) return(data.frame())
  fields <- names(rows[[1L]])
  out <- do.call(rbind, lapply(rows, function(row) {
    values <- lapply(fields, function(field) {
      value <- row[[field]]
      if (is.null(value)) NA else if (length(value) == 1L) value else {
        jsonlite::toJSON(value, auto_unbox = TRUE)
      }
    })
    names(values) <- fields
    as.data.frame(values, stringsAsFactors = FALSE, optional = TRUE)
  }))
  rownames(out) <- NULL
  out
}

#' Write a complete benchmark configuration template
#'
#' Most users can start from a preset and direct arguments. This function is
#' available when a complete custom configuration is needed.
#'
#' @param path Destination YAML path. Existing files are not replaced.
#' @return The normalized path written by the Python backend.
#' @export
#' @examples
#' \dontrun{mb_write_config_template("benchmark.yaml")}
mb_write_config_template <- function(path) {
  path <- .mb_normalize_path(path, must_work = FALSE)
  backend <- .mb_ensure_backend()
  .mb_python_path(backend$write_config_template(path))
}

#' Generate a truth-informed benchmark dataset
#'
#' Choose exactly one scenario source: a built-in `preset` or a complete
#' `config`. The configuration may be a YAML path or a named R list with the
#' same nested structure. Common changes have direct R arguments. Advanced
#' preset changes use a named `overrides` list, either nested like the YAML or
#' using dotted paths.
#'
#' Small runs return all benchmark tables as R data frames. Large runs require
#' `output_dir` and are streamed by the Python backend; the result then contains
#' paths, configuration, feature truth, and the memory estimate without copying
#' all sample-by-taxon tables into R.
#'
#' @param model An `mbbench_model` object, or `NULL` when `bank` is used.
#' @param bank An `mbbench_bank` object, or `NULL` when `model` is used.
#' @param preset Built-in benchmark preset name. Use [mb_list_presets()] to
#'   inspect the available choices.
#' @param config Path to a complete custom benchmark YAML, or a named R list.
#'   Cannot be combined with `preset` or `overrides`.
#' @param n_samples Optional preset sample-count override.
#' @param n_batches Optional preset batch-count override.
#' @param seed Optional preset seed override.
#' @param batch_size Optional preset generation chunk-size override.
#' @param overrides Named list of additional nested settings or dotted YAML
#'   paths and values.
#' @param output_dir Optional benchmark bundle directory.
#' @param overwrite Replace a nonempty output directory.
#' @param load Whether to return benchmark tables in R when they fit within
#'   `max_in_memory_mb`. If `FALSE`, `output_dir` is required and all sample
#'   tables remain on disk. Large requests are streamed even when `TRUE`.
#' @param max_in_memory_mb Maximum estimated memory for returning full tables.
#' @return An `mbbench_benchmark` list with the same table names as
#'   [mb_dataset()], plus `files`, `output_dir`, `loaded`, and
#'   `estimated_memory_mb`.
#' @export
#' @examples
#' \dontrun{
#' model <- mb_load_model(path = "my_model_bundle")
#' bench <- mb_benchmark(
#'   model,
#'   preset = "da-strict-common-moderate",
#'   n_samples = 120
#' )
#' bench$feature_truth
#' }
mb_benchmark <- function(model = NULL, bank = NULL,
                         preset = NULL, config = NULL,
                         n_samples = NULL, n_batches = NULL, seed = NULL,
                         batch_size = NULL, overrides = list(),
                         output_dir = NULL, overwrite = FALSE,
                         load = TRUE,
                         max_in_memory_mb = 512) {
  source <- .mb_select_source(model, bank)
  if (sum(c(!is.null(preset), !is.null(config))) != 1L) {
    .mb_abort("Provide exactly one of preset or config.")
  }
  if (length(overrides)) overrides <- .mb_named_list(overrides, "overrides")
  if (!is.null(config) && (length(overrides) || any(!vapply(
      list(n_samples, n_batches, seed, batch_size), is.null, logical(1))))) {
    .mb_abort("A custom config is complete and cannot be combined with preset overrides.")
  }
  if (!is.null(n_samples)) {
    n_samples <- .mb_scalar_number(n_samples, "n_samples", minimum = 1, integer = TRUE)
  }
  if (!is.null(n_batches)) {
    n_batches <- .mb_scalar_number(n_batches, "n_batches", minimum = 1, integer = TRUE)
  }
  if (!is.null(seed)) {
    seed <- .mb_scalar_number(seed, "seed", integer = TRUE)
  }
  if (!is.null(batch_size)) {
    batch_size <- .mb_scalar_number(batch_size, "batch_size", minimum = 1, integer = TRUE)
  }
  max_in_memory_mb <- .mb_scalar_number(max_in_memory_mb, "max_in_memory_mb", minimum = 0)
  if (max_in_memory_mb <= 0) .mb_abort("max_in_memory_mb must be positive.")
  if (!is.null(output_dir)) output_dir <- .mb_normalize_path(output_dir, must_work = FALSE)
  if (!is.logical(load) || length(load) != 1L || is.na(load)) {
    .mb_abort("load must be TRUE or FALSE.")
  }
  if (!load && is.null(output_dir)) {
    .mb_abort("load = FALSE requires output_dir for streamed generation.")
  }
  overwrite <- .mb_scalar_logical(overwrite, "overwrite")
  backend <- .mb_ensure_backend()
  if (!is.null(preset)) preset <- .mb_scalar_string(preset, "preset")
  if (!is.null(config)) {
    if (is.character(config)) {
      config <- .mb_normalize_path(config, must_work = TRUE)
    } else {
      config <- .mb_named_list(config, "config")
    }
  }
  arguments <- list(
    preset = preset, config = config, n_samples = n_samples,
    n_batches = n_batches, seed = seed, batch_size = batch_size,
    overrides = overrides, output_dir = output_dir,
    overwrite = overwrite, load = load,
    max_in_memory_mb = max_in_memory_mb
  )
  if (source$type == "model") arguments$model <- source$value$simulator
  else arguments$bank <- source$value$bank
  result <- do.call(backend$benchmark, arguments)
  .mb_benchmark_output(result)
}
