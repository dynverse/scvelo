#' Add the future dimensionality reduction based on RNA velocity information, using the existing dimensionality reduction
#'
#' @param dimred The dimensionality reduction
#' @export
add_dimred_future <- function(
  dataset,
  dimred = dynwrap::get_dimred(dataset),
  expression = dataset$expression,
  expression_future = dataset$expression_future
) {
  dataset$dimred_future <- embed_velocity(
    dataset,
    dimred = dimred,
    expression = expression,
    expression_future = expression_future
  )
  attr(dataset$dimred_future, "dimred_digest") <- digest::digest(dimred, algo = "md5")
  dataset
}

embed_velocity <- function(
  dataset,
  dimred = dynwrap::get_dimred(dataset),
  expression = dataset$expression,
  expression_future = dataset$expression_future
) {
  if(check_scvelo(dataset$velocity$scvelo)) {
    dimred_future <- embed_velocity_scvelo(
      dataset$velocity$scvelo,
      dimred = dimred
    )
  } else {
    dimred_future <- embed_velocity_new(
      dataset,
      dimred = dimred,
      expression = expression,
      expression_future = expression_future
    )
  }

  dimred_future
}


embed_velocity_new <- function(
  dataset,
  dimred = dynwrap::get_dimred(dataset),
  expression = dataset$expression,
  expression_future = dataset$expression_future
) {
  assertthat::assert_that(!is.null(expression))
  assertthat::assert_that(!is.null(expression_future))
  assertthat::assert_that(!is.null(dimred))

  # create adata object
  velo <- as.matrix(expression_future - expression)

  adata = anndata$AnnData(expression)
  adata$var_names <- colnames(expression)
  adata$obs_names <- rownames(expression)

  py_assign(adata$layers, "spliced", expression)
  py_assign(adata$layers, "velo", velo)

  # is necessary internally
  adata$uns[["velo_settings"]] <- list(mode = "deterministic")

  # embed velocity
  scvelo$tl$velocity_graph(adata, vkey = "velo")

  # assign dimred and embed
  py_assign(adata$obsm, "X_dimred", dimred)
  scvelo$tl$velocity_embedding(adata, basis = "dimred", vkey = "velo")

  dimred_future <- dimred + adata$obsm[["velo_dimred"]]
  dimnames(dimred_future) <- dimnames(dimred)

  # replace the ocasional NA dimred position
  dimred_future[is.na(dimred_future)] <- dimred[is.na(dimred_future)]

  dimred_future
}



# embed the velocity directly using an scvelo anndata object, instead of recreating one
embed_velocity_scvelo <- function(adata, dimred) {
  py_assign(adata$obsm, "X_dimred", dimred)
  scvelo$tl$velocity_embedding(adata, "dimred")
  dimred_future <- dimred + adata$obsm[["velocity_dimred"]]
  dimnames(dimred_future) <- dimnames(dimred)

  # replace the ocasional NA dimred position
  dimred_future[is.na(dimred_future)] <- dimred[is.na(dimred_future)]

  dimred_future
}
