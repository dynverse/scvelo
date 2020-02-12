#' Wrapper for the awesome scvelo package
#'
#' @import reticulate
#' @importFrom assertthat assert_that
#' @importFrom purrr map pmap_dbl
#'
#' @docType package
#' @name scelo
NULL

# global reference to scipy (will be initialized in .onLoad)
scvelo <- NULL
scanpy <- NULL
anndata <- NULL

.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scipy
  scvelo <<- import("scvelo", delay_load = TRUE)
  scanpy <<- import("scanpy", delay_load = TRUE)
  anndata <<- import("anndata", delay_load = TRUE)

  if(!py_module_available("scvelo")) {
    message("The scvelo python package is not installed, install it using `scvelo::install_scvelo()`")
  }
}


install_scvelo <- function(method = "auto", conda = "auto") {
  py_install("scvelo", method = method, conda = conda)
}
