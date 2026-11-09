.mb_model_name <- function(model) {
  model <- match.arg(model, c("ae", "vae", "ae-system"))
  model
}

.mb_model_default_dir <- function(pretrained, name) {
  as.character(pretrained$model_cache_dir(name))
}

.mb_model_url <- function(pretrained, name) {
  urls <- reticulate::py_to_r(pretrained$OFFICIAL_MODEL_URLS)
  value <- urls[[name]] %||% ""
  as.character(value)
}

.mb_prompt_directory <- function(default) {
  if (!interactive()) return(default)
  answer <- readline(sprintf("Install the model in [%s]: ", default))
  if (nzchar(trimws(answer))) .mb_normalize_path(trimws(answer), must_work = FALSE) else default
}

.mb_model_object <- function(simulator, path, device) {
  config <- reticulate::py_to_r(simulator$bundle$config)
  structure(list(
    simulator = simulator,
    path = .mb_normalize_path(path, must_work = TRUE),
    model_type = as.character(simulator$bundle$model_type),
    condition_fields = as.character(unlist(config$condition_fields, use.names = FALSE)),
    supported_output_modes = as.character(unlist(config$supported_output_modes, use.names = FALSE)),
    device = as.character(simulator$device)
  ), class = c("mbbench_model", "list"))
}

#' Download an official mbBench model
#'
#' If `install_dir` is omitted, an interactive session asks where to install the
#' model and displays the normal mbBench cache directory as the default.
#'
#' @param model Official model: `"ae"`, `"vae"`, or `"ae-system"`.
#' @param install_dir Optional installation directory.
#' @param url Optional ZIP or tar.gz URL for a mirror or custom archive.
#'   Defaults to the official Zenodo ZIP for the selected model.
#' @param device Device used when loading the returned model object.
#' @return An `mbbench_model` object ready for [mb_generate()] or
#'   [mb_benchmark()].
#' @export
#' @examples
#' \dontrun{model <- mb_download_model("ae")}
mb_download_model <- function(model = "ae",
                              install_dir = NULL, url = NULL, device = "auto") {
  name <- .mb_model_name(model)
  device <- .mb_scalar_string(device, "device")
  backend <- .mb_ensure_backend()
  pretrained <- reticulate::import("mbbench.pretrained", convert = FALSE)
  target <- install_dir %||% .mb_model_default_dir(pretrained, name)
  target <- .mb_normalize_path(target, must_work = FALSE)
  if (is.null(install_dir)) target <- .mb_prompt_directory(target)
  resolved_url <- if (is.null(url)) .mb_model_url(pretrained, name) else .mb_scalar_string(url, "url")
  if (!nzchar(resolved_url)) {
    .mb_abort(paste0(
      "No official download URL is registered for model '%s' in this release. ",
      "Provide url explicitly or load a model bundle with mb_load_model(path = PATH)."
    ), name)
  }
  pretrained$download_model(url = resolved_url, target_dir = target,
                            progress = .mb_show_download_progress())
  simulator <- backend$Simulator$from_directory(target, device = device)
  .mb_model_object(simulator, target, device)
}

#' Load a model bundle for repeated generation
#'
#' Supply `path` for a custom or manually downloaded model bundle. If `path` is
#' omitted, mbBench looks for the selected official model in its default cache;
#' an interactive session offers to download it when missing.
#'
#' The returned object keeps the Python simulator reference, so later calls only
#' need `model = model` and do not repeat the model directory.
#'
#' @param path Optional model bundle directory.
#' @param model Official model name used when `path` is omitted.
#' @param device `"auto"`, `"cpu"`, or a device understood by the Python backend.
#' @return An `mbbench_model` list.
#' @export
#' @examples
#' \dontrun{
#' model <- mb_load_model(path = "my_model_bundle")
#' model
#' }
mb_load_model <- function(path = NULL, model = "ae", device = "auto") {
  device <- .mb_scalar_string(device, "device")
  backend <- .mb_ensure_backend()
  if (!is.null(path)) {
    path <- .mb_normalize_path(path, must_work = TRUE)
    simulator <- backend$Simulator$from_directory(path, device = device)
    return(.mb_model_object(simulator, path, device))
  }
  name <- .mb_model_name(model)
  pretrained <- reticulate::import("mbbench.pretrained", convert = FALSE)
  target <- .mb_model_default_dir(pretrained, name)
  if (dir.exists(target)) {
    simulator <- backend$Simulator$from_directory(target, device = device)
    return(.mb_model_object(simulator, target, device))
  }
  message <- sprintf(
    "Official model '%s' is not installed in %s. Use mb_load_model(path = PATH) for a custom bundle.",
    name, target
  )
  if (!interactive()) .mb_abort("%s", message)
  answer <- utils::askYesNo(sprintf("%s Download it now?", message), default = FALSE)
  if (!isTRUE(answer)) .mb_abort("%s", message)
  mb_download_model(name, device = device)
}

.mb_check_model <- function(model) {
  if (!inherits(model, "mbbench_model") || is.null(model$simulator)) {
    .mb_abort("model must be an object returned by mb_load_model() or mb_download_model().")
  }
  model
}

#' List condition values supported by a loaded model or bank
#'
#' @param model An `mbbench_model` object, or `NULL` when `bank` is used.
#' @param bank An `mbbench_bank` object, or `NULL` when `model` is used.
#' @return A named list. Each element contains accepted values for one model
#'   condition field.
#' @export
#' @examples
#' \dontrun{
#' model <- mb_load_model(path = "my_model_bundle")
#' mb_list_conditions(model)
#' }
mb_list_conditions <- function(model = NULL, bank = NULL) {
  source <- .mb_select_source(model, bank)
  fields <- source$value$condition_fields
  values <- lapply(fields, function(field) {
    if (source$type == "model") {
      source$value$simulator$bundle$condition_values(field)
    } else {
      source$value$bank$condition_values(field)
    }
  })
  values <- lapply(values, function(value) as.character(reticulate::py_to_r(value)))
  stats::setNames(values, fields)
}
