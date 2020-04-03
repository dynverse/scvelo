#' Wrapper for the awesome scvelo package
#'
#' "scvelo" allows you to estimate velocity vectors for each cell and gene in a single-cell expression dataset.
#'
#' @import reticulate
#' @import assertthat
#' @importFrom purrr %>% map pmap_dbl
#' @importFrom dplyr filter mutate arrange
#'
#' @docType package
#' @name scvelo
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

.onAttach <- function(libname, pkgname) {
  if (!py_module_available("scvelo")) {
    packageStartupMessage("The scvelo python package is not installed. Install it using `scvelo::install_scvelo()`")
  } else {
    vers <- package_version(gsub("\\.dev.*", "", scvelo$`__version__`))
    if (vers < "0.1.26") {
      packageStartupMessage("The scvelo python package should be >= 0.1.26. Install the scvelo from Github using `scvelo::install_scvelo()`")
    }
  }
}

install_scvelo <- function(method = "auto", conda = "auto") {
  message("Installing scvelo through bioconda currently does not work, will install scvelo from Github using pip.")
  # py_install("scvelo", method = "pip")
  reticulate::py_install(packages = "git+https://github.com/theislab/scvelo.git", pip = TRUE)
}
