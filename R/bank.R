.mb_bank_object <- function(bank, path) {
  manifest <- reticulate::py_to_r(bank$manifest)
  shape <- reticulate::py_to_r(bank$metadata$shape)
  structure(list(
    bank = bank,
    path = .mb_normalize_path(path, must_work = TRUE),
    bank_type = as.character(manifest$bank_type),
    condition_fields = as.character(reticulate::py_to_r(bank$condition_fields)),
    n_samples = as.integer(shape[[1L]]),
    n_taxa = length(reticulate::py_to_r(bank$taxa))
  ), class = c("mbbench_bank", "list"))
}

#' Load a batch-neutral logit bank
#'
#' A bank contains pre-generated abundance and support logits. Generation and
#' benchmark construction then run on CPU without sampling a diffusion model.
#'
#' @param path Directory containing `bank.yaml`, both logit matrices, metadata,
#'   and taxa.
#' @return An `mbbench_bank` object for [mb_generate()] or [mb_benchmark()].
#' @export
#' @examples
#' \dontrun{
#' bank <- mb_load_bank("my_bank")
#' generated <- mb_generate(bank = bank, n_samples = 20)
#' }
mb_load_bank <- function(path) {
  path <- .mb_normalize_path(path, must_work = TRUE)
  backend <- .mb_ensure_backend()
  .mb_bank_object(backend$load_bank(path), path)
}

#' Download and load a pre-generated logit bank
#'
#' @param name Name of the official bank. Currently only `core` is supported.
#' @param install_dir Optional installation directory. An existing valid bank
#'   at that path is reused.
#' @param url Optional direct ZIP URL, including a local `file://` URL.
#'   Defaults to the official core bank ZIP on Zenodo.
#' @return An `mbbench_bank` object for [mb_generate()] or [mb_benchmark()].
#' @export
mb_download_bank <- function(name = "core", install_dir = NULL, url = NULL) {
  name <- .mb_scalar_string(name, "name")
  if (!is.null(install_dir)) {
    install_dir <- .mb_normalize_path(.mb_scalar_string(install_dir, "install_dir"))
  }
  if (!is.null(url)) url <- .mb_scalar_string(url, "url")
  backend <- .mb_ensure_backend()
  bank <- backend$download_bank(name = name, install_dir = install_dir, url = url,
                                progress = .mb_show_download_progress())
  .mb_bank_object(bank, .mb_python_path(bank$root))
}

.mb_check_bank <- function(bank) {
  if (!inherits(bank, "mbbench_bank") || is.null(bank$bank)) {
    .mb_abort("bank must be an object returned by mb_load_bank().")
  }
  bank
}

.mb_select_source <- function(model, bank) {
  if (sum(c(!is.null(model), !is.null(bank))) != 1L) {
    .mb_abort("Provide exactly one of model or bank.")
  }
  if (!is.null(model)) {
    return(list(type = "model", value = .mb_check_model(model)))
  }
  list(type = "bank", value = .mb_check_bank(bank))
}
