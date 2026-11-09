.mb_read_json <- function(path) {
  value <- jsonlite::read_json(path, simplifyVector = FALSE)
  if (!is.list(value) || is.null(names(value))) {
    .mb_abort("JSON root must be an object: %s", path)
  }
  value
}

.mb_write_json <- function(path, value) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(value, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
}

.mb_header <- function(headers, name) {
  if (is.null(headers) || !length(headers)) return(NULL)
  keys <- tolower(names(headers))
  hit <- which(keys == tolower(name))
  if (!length(hit)) NULL else as.character(headers[[hit[[1L]]]])
}

.mb_cached_json <- function(url, name) {
  cache <- .mb_cache_dir("catalog")
  payload_path <- file.path(cache, paste0(name, ".json"))
  state_path <- file.path(cache, paste0(name, ".http.json"))
  state <- if (file.exists(state_path)) .mb_read_json(state_path) else list()
  handle <- curl::new_handle()
  curl::handle_setopt(handle, timeout = 30, followlocation = TRUE, failonerror = FALSE,
                     accept_encoding = "gzip")
  request_headers <- character()
  if (identical(state$url, url)) {
    if (!is.null(state$etag) && nzchar(state$etag)) {
      request_headers["If-None-Match"] <- state$etag
    }
    if (!is.null(state$last_modified) && nzchar(state$last_modified)) {
      request_headers["If-Modified-Since"] <- state$last_modified
    }
  }
  if (length(request_headers)) {
    do.call(curl::handle_setheaders, c(list(handle = handle), as.list(request_headers)))
  }

  response <- tryCatch(
    curl::curl_fetch_memory(url, handle = handle),
    error = identity
  )
  if (inherits(response, "error")) {
    if (file.exists(payload_path)) {
      warning(sprintf("Could not refresh the mbBench dataset catalog (%s); using the cached copy.",
                      conditionMessage(response)), call. = FALSE)
      return(.mb_read_json(payload_path))
    }
    .mb_abort(paste0(
      "Could not reach the mbBench dataset catalog. Download the dataset ZIP ",
      "manually and use mb_dataset(archive = PATH)."
    ))
  }
  if (identical(response$status_code, 304L) && file.exists(payload_path)) {
    return(.mb_read_json(payload_path))
  }
  if (response$status_code < 200L || response$status_code >= 300L) {
    if (file.exists(payload_path)) {
      warning(sprintf("Could not refresh the mbBench dataset catalog (HTTP %s); using the cached copy.",
                      response$status_code), call. = FALSE)
      return(.mb_read_json(payload_path))
    }
    .mb_abort("Could not download mbBench dataset metadata (HTTP %s).", response$status_code)
  }

  text <- rawToChar(response$content)
  payload <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  if (!is.list(payload) || is.null(names(payload))) {
    .mb_abort("Remote JSON root must be an object: %s", url)
  }
  dir.create(cache, recursive = TRUE, showWarnings = FALSE)
  writeLines(text, payload_path, useBytes = TRUE)
  headers <- tryCatch(curl::parse_headers_list(response$headers), error = function(e) list())
  .mb_write_json(state_path, list(
    url = url,
    etag = .mb_header(headers, "etag"),
    last_modified = .mb_header(headers, "last-modified")
  ))
  payload
}

.mb_catalog_url <- function(record) {
  files <- record$files
  entries <- if (is.list(files) && !is.null(files$entries)) files$entries else files
  for (entry in entries %||% list()) {
    if (!is.list(entry)) next
    filename <- entry$key %||% entry$filename
    if (!identical(filename, .mb_catalog_filename)) next
    links <- entry$links %||% list()
    url <- links$content %||% links$download %||% links$self
    if (!is.null(url) && nzchar(url)) return(url)
  }
  .mb_abort("The official dataset record does not contain %s.", .mb_catalog_filename)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.mb_validate_catalog <- function(catalog) {
  if (!identical(catalog$schema_version, .mb_catalog_schema)) {
    .mb_abort("Unsupported dataset catalog schema: %s.", catalog$schema_version %||% "<missing>")
  }
  if (!is.character(catalog$latest) || length(catalog$latest) != 1L || !nzchar(catalog$latest)) {
    .mb_abort("Dataset catalog requires a nonempty latest release.")
  }
  releases <- catalog$releases
  if (!is.list(releases) || !length(releases)) {
    .mb_abort("Dataset catalog requires at least one release.")
  }
  versions <- character()
  for (release in releases) {
    version <- release$version
    datasets <- release$datasets
    if (!is.character(version) || length(version) != 1L || !nzchar(version) || version %in% versions) {
      .mb_abort("Dataset release versions must be unique nonempty strings.")
    }
    versions <- c(versions, version)
    if (!is.list(datasets) || !length(datasets)) {
      .mb_abort("Dataset release '%s' contains no datasets.", version)
    }
    identities <- character()
    for (entry in datasets) {
      preset <- entry$preset
      seed <- entry$seed
      archive <- entry$archive
      url <- entry$url
      if (!is.character(preset) || length(preset) != 1L || !nzchar(preset)) {
        .mb_abort("Every dataset entry requires a preset name.")
      }
      if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed != as.integer(seed)) {
        .mb_abort("Dataset '%s' requires an integer seed.", preset)
      }
      if (!is.character(archive) || length(archive) != 1L || !grepl("\\.zip$", archive, ignore.case = TRUE)) {
        .mb_abort("Dataset '%s' requires a ZIP archive name.", preset)
      }
      if (!is.character(url) || length(url) != 1L || !nzchar(url)) {
        .mb_abort("Dataset '%s' requires a download URL.", preset)
      }
      identity <- paste(preset, as.integer(seed), sep = "\r")
      if (identity %in% identities) {
        .mb_abort("Dataset release '%s' repeats preset '%s', seed %s.", version, preset, seed)
      }
      identities <- c(identities, identity)
    }
  }
  if (!catalog$latest %in% versions) {
    .mb_abort("Latest dataset release '%s' is not listed.", catalog$latest)
  }
  invisible(catalog)
}

.mb_load_catalog <- function() {
  record_url <- getOption("mbBench.dataset_record_url", .mb_official_dataset_record_url)
  if (!is.character(record_url) || length(record_url) != 1L || !nzchar(record_url)) {
    .mb_abort(paste0(
      "The official mbBench dataset record has not been published yet. ",
      "Use mb_dataset(archive = PATH) for a manually downloaded ZIP."
    ))
  }
  record <- .mb_cached_json(record_url, "record")
  latest_url <- record$links$latest %||% NULL
  if (!is.null(latest_url) && nzchar(latest_url) && !identical(latest_url, record_url)) {
    record <- .mb_cached_json(latest_url, "latest_record")
  }
  catalog <- .mb_cached_json(.mb_catalog_url(record), "datasets")
  .mb_validate_catalog(catalog)
  catalog
}

.mb_release <- function(catalog, version) {
  selected <- if (identical(version, "latest")) catalog$latest else version
  for (release in catalog$releases) {
    if (identical(release$version, selected)) return(release)
  }
  available <- vapply(catalog$releases, `[[`, character(1), "version")
  .mb_abort("Unknown dataset version '%s'. Available versions: %s.",
            version, paste(available, collapse = ", "))
}

.mb_dataset_entries <- function(version = "latest") {
  version <- .mb_scalar_string(version, "version")
  release <- .mb_release(.mb_load_catalog(), version)
  lapply(release$datasets, function(entry) {
    c(entry, list(
      version = release$version,
      version_doi = entry$artifact_version_doi %||% release$version_doi,
      record_id = entry$artifact_record_id %||% release$record_id
    ))
  })
}

#' List pre-generated benchmark datasets
#'
#' Reads the official catalog and returns one row per preset and seed. The
#' catalog is checked automatically when this function is called.
#'
#' @param version Dataset release version, or `"latest"`.
#' @return A data frame with catalog fields such as `preset`, `seed`, `version`,
#'   `n_samples`, and `n_batches` when supplied by the release.
#' @export
#' @examples
#' \dontrun{mb_list_datasets()}
mb_list_datasets <- function(version = "latest") {
  entries <- .mb_dataset_entries(version)
  fields <- unique(unlist(lapply(entries, names), use.names = FALSE))
  rows <- lapply(entries, function(entry) {
    values <- stats::setNames(vector("list", length(fields)), fields)
    for (field in fields) values[[field]] <- entry[[field]] %||% NA
    as.data.frame(values, stringsAsFactors = FALSE, optional = TRUE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Show one pre-generated benchmark dataset
#'
#' @param preset Preset name.
#' @param version Dataset release version, or `"latest"`.
#' @return A data frame with one row per published seed.
#' @export
#' @examples
#' \dontrun{mb_show_dataset("da-strict-common-moderate")}
mb_show_dataset <- function(preset, version = "latest") {
  preset <- .mb_scalar_string(preset, "preset")
  result <- mb_list_datasets(version)
  matches <- result[result$preset == preset, , drop = FALSE]
  if (!nrow(matches)) {
    available <- unique(result$preset)
    .mb_abort("Unknown pre-generated dataset '%s'. Available presets: %s.",
              preset, paste(available, collapse = ", "))
  }
  rownames(matches) <- NULL
  matches
}

.mb_read_csv_gz <- function(path) {
  connection <- gzfile(path, open = "rt", encoding = "UTF-8")
  on.exit(close(connection), add = TRUE)
  utils::read.csv(connection, check.names = FALSE, stringsAsFactors = FALSE)
}

.mb_bundle_files <- c(
  counts = "counts.csv.gz",
  relative_abundance = "relative_abundance.csv.gz",
  latent_composition = "latent_composition.csv.gz",
  metadata = "metadata.csv.gz",
  feature_truth = "feature_truth.csv.gz",
  sample_truth = "sample_truth.csv.gz",
  batch_neutral_reference = "batch_neutral_biological_reference.csv.gz"
)

.mb_validate_bundle <- function(tables, config) {
  taxa <- setdiff(names(tables$counts), "sample_id")
  if (!length(taxa) || !"sample_id" %in% names(tables$counts)) {
    .mb_abort("Benchmark counts require a sample_id column and at least one taxon.")
  }
  if (any(!nzchar(trimws(taxa))) || anyDuplicated(taxa)) {
    .mb_abort("Benchmark taxon names must be nonempty and unique.")
  }
  abundance_names <- c("relative_abundance", "latent_composition", "batch_neutral_reference")
  for (name in abundance_names) {
    if (!identical(names(tables[[name]]), c("sample_id", taxa))) {
      .mb_abort("Benchmark abundance tables must have the same ordered taxa.")
    }
  }
  sample_ids <- as.character(tables$counts$sample_id)
  if (anyNA(tables$counts$sample_id) || any(!nzchar(trimws(sample_ids))) ||
      anyDuplicated(sample_ids)) {
    .mb_abort("Benchmark sample IDs must be nonempty and unique.")
  }
  for (name in c(abundance_names, "metadata", "sample_truth")) {
    if (!"sample_id" %in% names(tables[[name]]) ||
        !identical(as.character(tables[[name]]$sample_id), sample_ids)) {
      .mb_abort("%s sample IDs or row order do not match counts.", name)
    }
  }
  counts <- as.matrix(tables$counts[taxa])
  storage.mode(counts) <- "double"
  depths <- rowSums(counts)
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts)) || any(depths <= 0)) {
    .mb_abort("Benchmark counts must be nonnegative integers with positive row depth.")
  }
  relative <- as.matrix(tables$relative_abundance[taxa])
  if (!isTRUE(all.equal(relative, counts / depths, tolerance = 1e-12, check.attributes = FALSE))) {
    .mb_abort("relative_abundance must equal counts divided by row depth.")
  }
  for (name in c("latent_composition", "batch_neutral_reference")) {
    values <- as.matrix(tables[[name]][taxa])
    if (any(!is.finite(values)) || any(values < 0) || any(abs(rowSums(values) - 1) > 1e-10)) {
      .mb_abort("%s must be finite, nonnegative, and sum to one by row.", name)
    }
  }
  if (!"taxon" %in% names(tables$feature_truth) ||
      !identical(as.character(tables$feature_truth$taxon), taxa)) {
    .mb_abort("feature_truth taxon order does not match abundance tables.")
  }
  if (!is.null(tables$batch_neutral_counts)) {
    if (!identical(names(tables$batch_neutral_counts), c("sample_id", taxa)) ||
        !identical(as.character(tables$batch_neutral_counts$sample_id), sample_ids)) {
      .mb_abort("batch_neutral_counts does not match the primary tables.")
    }
    neutral <- as.matrix(tables$batch_neutral_counts[taxa])
    storage.mode(neutral) <- "double"
    if (any(!is.finite(neutral)) || any(neutral < 0) || any(neutral != round(neutral))) {
      .mb_abort("batch_neutral_counts must contain nonnegative integers.")
    }
  }
  invisible(config)
}

.mb_load_bundle <- function(path, source = NULL) {
  root <- .mb_normalize_path(path, must_work = TRUE)
  if (!dir.exists(root)) .mb_abort("Dataset directory does not exist: %s", root)
  expected <- c("benchmark.yaml", unname(.mb_bundle_files))
  missing <- expected[!file.exists(file.path(root, expected))]
  if (length(missing)) {
    .mb_abort("Benchmark dataset is incomplete; missing: %s.", paste(missing, collapse = ", "))
  }
  tables <- lapply(.mb_bundle_files, function(filename) .mb_read_csv_gz(file.path(root, filename)))
  neutral <- file.path(root, "batch_neutral_counts.csv.gz")
  tables$batch_neutral_counts <- if (file.exists(neutral)) .mb_read_csv_gz(neutral) else NULL
  config <- yaml::read_yaml(file.path(root, "benchmark.yaml"))
  .mb_validate_bundle(tables, config)
  structure(c(tables, list(config = config, path = root, source = source)),
            class = c("mbbench_dataset", "list"))
}

.mb_find_bundle_root <- function(path) {
  configs <- list.files(path, pattern = "^benchmark\\.yaml$", recursive = TRUE,
                        full.names = TRUE, include.dirs = FALSE)
  if (length(configs) != 1L) {
    .mb_abort("Dataset archive must contain exactly one mbBench benchmark bundle.")
  }
  dirname(configs[[1L]])
}

.mb_install_archive <- function(archive, target, overwrite) {
  archive <- .mb_normalize_path(archive, must_work = TRUE)
  target <- .mb_normalize_path(target, must_work = FALSE)
  if (dir.exists(target) && !overwrite) return(target)
  extracted <- tempfile("mbbench-unzip-")
  dir.create(extracted)
  on.exit(unlink(extracted, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(archive, exdir = extracted)
  source <- .mb_find_bundle_root(extracted)
  .mb_load_bundle(source)
  if (dir.exists(target)) unlink(target, recursive = TRUE, force = TRUE)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  files <- list.files(source, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!all(file.copy(files, target, recursive = TRUE, copy.date = TRUE))) {
    unlink(target, recursive = TRUE, force = TRUE)
    .mb_abort("Could not install dataset archive in %s.", target)
  }
  target
}

.mb_import_dataset <- function(archive, install_dir = NULL, overwrite = FALSE) {
  archive <- .mb_normalize_path(archive, must_work = TRUE)
  overwrite <- .mb_scalar_logical(overwrite, "overwrite")
  if (!file.exists(archive) || dir.exists(archive)) {
    .mb_abort("Dataset archive does not exist: %s", archive)
  }
  target <- install_dir %||% file.path(.mb_cache_dir("datasets", "manual"),
                                       tools::file_path_sans_ext(basename(archive)))
  root <- .mb_install_archive(archive, target, overwrite)
  .mb_load_bundle(root, source = list(type = "manual_archive", archive = archive))
}

#' Download and load a pre-generated benchmark dataset
#'
#' The selected ZIP is installed in the mbBench user cache by default and all
#' tables are then loaded into R. Set `install_dir` to keep it elsewhere.
#'
#' @param preset Preset name from [mb_list_datasets()].
#' @param seed Published generation seed.
#' @param version Dataset release version, or `"latest"`.
#' @param install_dir Optional installation directory.
#' @param overwrite Replace an existing installation.
#' @return An `mbbench_dataset` list. See [mb_dataset()].
#' @export
#' @examples
#' \dontrun{dat <- mb_download_dataset("da-strict-common-moderate")}
mb_download_dataset <- function(preset, seed = 42L, version = "latest",
                                install_dir = NULL, overwrite = FALSE) {
  preset <- .mb_scalar_string(preset, "preset")
  seed <- .mb_scalar_number(seed, "seed", integer = TRUE)
  overwrite <- .mb_scalar_logical(overwrite, "overwrite")
  entries <- .mb_dataset_entries(version)
  matches <- Filter(function(x) identical(x$preset, preset) && as.integer(x$seed) == seed, entries)
  if (length(matches) != 1L) {
    choices <- vapply(entries, function(x) sprintf("%s (seed %s)", x$preset, x$seed), character(1))
    .mb_abort("No dataset matches preset '%s', seed %s. Available datasets: %s.",
              preset, seed, paste(choices, collapse = ", "))
  }
  entry <- matches[[1L]]
  target <- install_dir %||% file.path(.mb_cache_dir("datasets"), entry$version,
                                       sprintf("%s_seed%s", preset, seed))
  if (dir.exists(target) && !overwrite) {
    root <- .mb_normalize_path(target, must_work = TRUE)
  } else {
    archive <- tempfile(fileext = ".zip")
    on.exit(unlink(archive, force = TRUE), add = TRUE)
    error <- tryCatch({
      curl::curl_download(entry$url, archive, quiet = !.mb_show_download_progress())
      NULL
    }, error = identity)
    if (inherits(error, "error")) {
      .mb_abort(paste0(
        "Could not download %s. Download the ZIP manually and use ",
        "mb_dataset(archive = PATH)."
      ), entry$archive)
    }
    root <- .mb_install_archive(archive, target, overwrite = TRUE)
  }
  .mb_load_bundle(root, source = list(
    type = "zenodo_download",
    requested_version = version,
    resolved_version = entry$version,
    version_doi = entry$version_doi %||% NULL,
    record_id = entry$record_id %||% NULL,
    preset = preset,
    seed = seed
  ))
}

#' Load a pre-generated mbBench dataset
#'
#' Provide exactly one source: an official `preset`, a manually downloaded ZIP
#' `archive`, or an extracted bundle `path`. Every table is fully loaded into R;
#' this function does not install or initialize the Python backend.
#'
#' The returned list contains `counts`, `relative_abundance`,
#' `latent_composition`, `metadata`, `feature_truth`, `sample_truth`,
#' `batch_neutral_reference`, optional `batch_neutral_counts`, `config`, `path`,
#' and `source`.
#'
#' @param preset Optional official preset name.
#' @param archive Optional path to a manually downloaded ZIP.
#' @param path Optional path to an extracted benchmark bundle.
#' @param seed Published generation seed used with `preset`.
#' @param version Dataset release version used with `preset`.
#' @param install_dir Optional installation directory for a ZIP or download.
#' @param overwrite Replace an existing installation.
#' @return An `mbbench_dataset` object, which is an ordinary named list of R
#'   data frames and configuration metadata.
#' @export
#' @examples
#' \dontrun{
#' dat <- mb_dataset(path = "benchmark_bundle")
#' dat$counts
#' dat$feature_truth
#' }
mb_dataset <- function(preset = NULL, archive = NULL, path = NULL, seed = 42L,
                       version = "latest", install_dir = NULL, overwrite = FALSE) {
  overwrite <- .mb_scalar_logical(overwrite, "overwrite")
  supplied <- c(!is.null(preset), !is.null(archive), !is.null(path))
  if (sum(supplied) != 1L) {
    .mb_abort("Provide exactly one of preset, archive, or path.")
  }
  if (!is.null(preset)) {
    return(mb_download_dataset(preset, seed, version, install_dir, overwrite))
  }
  if (!is.null(archive)) {
    return(.mb_import_dataset(archive, install_dir, overwrite))
  }
  .mb_load_bundle(path, source = list(type = "local_directory",
                                      path = .mb_normalize_path(path, must_work = TRUE)))
}
