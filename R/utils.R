.mb_abort <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

.mb_scalar_string <- function(x, name, allow_empty = FALSE) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      (!allow_empty && !nzchar(x))) {
    .mb_abort("%s must be one non-missing string.", name)
  }
  x
}

.mb_scalar_number <- function(x, name, minimum = -Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < minimum || (integer && x != as.integer(x))) {
    qualifier <- if (integer) "integer" else "number"
    .mb_abort("%s must be one finite %s not smaller than %s.", name, qualifier, minimum)
  }
  if (integer) as.integer(x) else as.numeric(x)
}

.mb_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .mb_abort("%s must be TRUE or FALSE.", name)
  }
  isTRUE(x)
}

.mb_cache_dir <- function(...) {
  override <- getOption("mbBench.cache_dir", NULL)
  root <- if (is.null(override)) tools::R_user_dir("mbBench", "cache") else override
  file.path(root, ...)
}

.mb_normalize_path <- function(path, must_work = FALSE) {
  normalizePath(path.expand(path), winslash = "/", mustWork = must_work)
}

.mb_show_download_progress <- function() {
  interactive() || isatty(stderr())
}

.mb_named_list <- function(x, name) {
  if (!is.list(x) || is.null(names(x)) || any(!nzchar(names(x)))) {
    .mb_abort("%s must be a named list.", name)
  }
  duplicates <- unique(names(x)[duplicated(names(x))])
  if (length(duplicates)) {
    .mb_abort(
      "%s contains duplicate names: %s.",
      name,
      paste(duplicates, collapse = ", ")
    )
  }
  x
}

.mb_require_namespace <- function(package, reason) {
  if (!requireNamespace(package, quietly = TRUE)) {
    .mb_abort("Package '%s' is required %s. Install it with install.packages(\"%s\").",
              package, reason, package)
  }
  invisible(TRUE)
}
