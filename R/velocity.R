#' Add velocity to a dynwrap dataset
#'
#' @inheritParams dynwrap::add_trajectory
#'
#' @export
add_velocity <- function(
  dataset,
  spliced = dataset$expression,
  unspliced = dataset$expression_unspliced,
  mode = "stochastic",
  velocity = get_velocity(spliced, unspliced, mode)
) {
  assertthat::assert_that(!is.null(velocity$expression_projected))

  dataset$expression_projected <- velocity$expression_projected

  velocity$expression_projected <- NULL
  dataset$velocity <- velocity

  dynutils::add_class(dataset, "wrapper_with_velocity")

  dataset
}


#' Calculate velocity
#'
#' @param spliced Spliced expression matrix
#' @param unspliced Unspliced expression matrix
#'
#' @export
get_velocity <- function(spliced, unspliced, mode = "stochastic") {
  assertthat::assert_that(all(dim(spliced) == dim(unspliced)))
  assertthat::assert_that(all(rownames(spliced) == rownames(unspliced)))
  assertthat::assert_that(all(colnames(spliced) == colnames(unspliced)))

  # create anndata object
  velocity = anndata$AnnData(spliced)
  velocity$var_names <- colnames(spliced)
  velocity$obs_names <- rownames(spliced)

  py_assign(velocity$layers, "spliced", spliced)
  py_assign(velocity$layers, "unspliced", unspliced)

  # calculate velocity
  # py_capture_output({ # can't capture output because of https://github.com/rstudio/reticulate/issues/386, otherwise crash when testing

    scvelo$pp$moments(velocity, n_neighbors = 20L)
    if (mode == "dynamical") {
      scvelo$tl$recover_dynamics(velocity)
    }
    scvelo$tl$velocity(velocity, mode=mode)
    scvelo$tl$velocity_graph(velocity)
  # })

  # return vecloity object
  tibble::lst(
    expression_projected = spliced + velocity$layers[["velocity"]],
    scvelo = velocity
  )
}

#' Embed the velocity within a given dimensionality reduction
#'
#' @param dimred The dimensionality reduction
#' @export
embed_velocity <- function(
  dataset,
  expression = dataset$expression,
  expression_projected = dataset$expression_projected,
  dimred = dynwrap::get_dimred(dataset)
) {
  assertthat::assert_that(!is.null(expression))
  assertthat::assert_that(!is.null(expression_projected))
  assertthat::assert_that(!is.null(dimred))

  # create adata object
  velo <- as.matrix(expression_projected - expression)

  adata = anndata$AnnData(expression)
  adata$var_names <- colnames(expression)
  adata$obs_names <- rownames(expression)

  py_assign(adata$layers, "spliced", expression)
  py_assign(adata$layers, "velo", velo)

  # is necessary internally
  adata$uns[["velo_settings"]] <- list(mode = "deterministic")

  # assign dimred
  py_assign(adata$obsm, "X_dimred", as.matrix(dimred))

  # embed velocity
  scvelo$tl$velocity_graph(adata, vkey = "velo")
  scvelo$tl$velocity_embedding(adata, vkey = "velo", basis = "dimred")

  # assign dimred and embed
  py_assign(adata$obsm, "X_dimred", dimred)
  scvelo$tl$velocity_embedding(adata, basis = "dimred", vkey = "velo")

  dimred_projected <- dimred + adata$obsm[["velo_dimred"]]
  dimnames(dimred_projected) <- dimnames(dimred)

  # replace the ocasional NA dimred position
  dimred_projected[is.na(dimred_projected)] <- dimred[is.na(dimred_projected)]

  dimred_projected
}



# embed the velocity directly using an scvelo anndata object, instead of recreating one
embed_velocity_scvelo <- function(adata, dimred) {
  py_assign(adata$obsm, "X_dimred", dimred)
  scvelo$tl$velocity_embedding(adata, "dimred")
  dimred_projected <- dimred + adata$obsm[["velocity_dimred"]]
  dimnames(dimred_projected) <- dimnames(dimred)

  # replace the ocasional NA dimred position
  dimred_projected[is.na(dimred_projected)] <- dimred[is.na(dimred_projected)]

  dimred_projected
}
