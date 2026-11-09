.mb_backend_config_path <- function() .mb_cache_dir("backend.json")

.mb_read_backend_config <- function() {
  path <- .mb_backend_config_path()
  if (!file.exists(path)) return(NULL)
  .mb_read_json(path)
}

.mb_backend_spec <- function(source = NULL) {
  if (is.null(source)) return(.mb_backend_requirement)
  source <- .mb_scalar_string(source, "source")
  if (file.exists(path.expand(source))) {
    path <- .mb_normalize_path(source, must_work = TRUE)
    if (!dir.exists(path)) .mb_abort("Local backend source must be a package directory: %s", path)
    return(paste0("mbBench @ file:///", utils::URLencode(path)))
  }
  source
}

.mb_declare_backend <- function(spec, python_version) {
  .mb_require_namespace("reticulate", "to use the Python generation backend")
  if (utils::packageVersion("reticulate") < "1.41") {
    .mb_abort("reticulate 1.41 or newer is required. Update it with install.packages(\"reticulate\").")
  }
  reticulate::py_require(packages = spec, python_version = python_version)
  invisible(TRUE)
}

.mb_import_backend <- function() {
  reticulate::import("mbbench", convert = FALSE)
}

#' Install or configure the mbBench Python backend
#'
#' Declares the Python dependency with `reticulate`, lets `reticulate` construct
#' its managed environment, and imports the backend once to verify it. Call this
#' function explicitly, or accept the prompt shown by the first generation call.
#'
#' Set `source` to a local Python package directory to use a source checkout.
#'
#' @param source `NULL` for the compatible PyPI release, a local Python package
#'   directory, or a pip-compatible requirement string.
#' @param python_version Python version constraint passed to
#'   [reticulate::py_require()].
#' @return An invisible backend status list.
#' @export
#' @examples
#' \dontrun{
#' mb_install_backend()
#' mb_install_backend(source = "../mbBench")
#' }
mb_install_backend <- function(source = NULL, python_version = ">=3.10,<3.13") {
  python_version <- .mb_scalar_string(python_version, "python_version")
  spec <- .mb_backend_spec(source)
  .mb_declare_backend(spec, python_version)
  module <- tryCatch(.mb_import_backend(), error = identity)
  if (inherits(module, "error")) {
    .mb_abort("Could not install or import the mbBench Python backend: %s",
              conditionMessage(module))
  }
  config <- list(requirement = spec, python_version = python_version)
  .mb_write_json(.mb_backend_config_path(), config)
  invisible(mb_backend_status(initialize = TRUE))
}

#' Report Python backend status
#'
#' By default this reads local configuration without starting Python. Use
#' `initialize = TRUE` to import and verify the configured backend.
#'
#' @param initialize Whether to initialize Python and verify the import.
#' @return A named list with `configured`, `initialized`, `checked`, `available`,
#'   `requirement`, `python`, and `version`. `available` is `NA` until the
#'   backend has actually been checked.
#' @export
#' @examples
#' mb_backend_status()
mb_backend_status <- function(initialize = FALSE) {
  config <- .mb_read_backend_config()
  has_reticulate <- requireNamespace("reticulate", quietly = TRUE)
  initialized <- has_reticulate && reticulate::py_available(initialize = FALSE)
  checked <- FALSE
  available <- NA
  python <- NULL
  version <- NULL
  if (isTRUE(initialize)) {
    .mb_require_namespace("reticulate", "to inspect the Python generation backend")
    if (!is.null(config)) {
      .mb_declare_backend(config$requirement, config$python_version)
    }
    module <- tryCatch(.mb_import_backend(), error = identity)
    available <- !inherits(module, "error")
    checked <- TRUE
    initialized <- reticulate::py_available(initialize = FALSE)
    if (available) version <- tryCatch(as.character(module$`__version__`), error = function(e) NULL)
  } else if (initialized) {
    available <- isTRUE(tryCatch(reticulate::py_module_available("mbbench"), error = function(e) FALSE))
    checked <- TRUE
    if (available) {
      module <- reticulate::import("mbbench", convert = FALSE)
      version <- tryCatch(as.character(module$`__version__`), error = function(e) NULL)
    }
  }
  if (initialized) {
    python <- tryCatch(reticulate::py_config()$python, error = function(e) NULL)
  }
  structure(list(
    configured = !is.null(config),
    initialized = initialized,
    checked = checked,
    available = available,
    requirement = config$requirement %||% .mb_backend_requirement,
    python = python,
    version = version
  ), class = c("mbbench_backend_status", "list"))
}

.mb_ensure_backend <- function() {
  .mb_require_namespace("reticulate", "to generate new data")
  config <- .mb_read_backend_config()
  if (!is.null(config)) {
    .mb_declare_backend(config$requirement, config$python_version)
    module <- tryCatch(.mb_import_backend(), error = identity)
    if (!inherits(module, "error")) return(module)
    .mb_abort(paste0(
      "The configured mbBench Python backend could not be imported: %s. ",
      "Restart R and run mb_install_backend() again, or select a working Python ",
      "before the first reticulate call."
    ), conditionMessage(module))
  }

  existing <- tryCatch(reticulate::py_module_available("mbbench"), error = function(e) FALSE)
  if (isTRUE(existing)) return(.mb_import_backend())
  if (!interactive()) {
    .mb_abort(paste0(
      "The mbBench Python backend is not configured. Run mb_install_backend() ",
      "before generating data. Loading pre-generated data with mb_dataset() does not need it."
    ))
  }
  answer <- utils::askYesNo("The mbBench Python backend is not installed. Install it now?",
                            default = FALSE)
  if (!isTRUE(answer)) {
    .mb_abort("Python backend installation was not selected. Run mb_install_backend() when ready.")
  }
  mb_install_backend()
  .mb_import_backend()
}
