#' Wrapper for the awesome scvelo package
#'
#' @import reticulate
#'
#' @docType package
#' @name reticulate
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
}


install_scvelo <- function(method = "auto", conda = "auto") {
  py_install("scvelo", method = method, conda = conda)
}
