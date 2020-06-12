#' Add the future dimensionality reduction based on RNA velocity information, using the existing dimensionality reduction
#'
#' @inheritParams dynwrap::add_expression
#' @inheritParams dynwrap::add_dimred
#'
#' @importFrom digest digest
#'
#' @export
add_dimred_future <- function(
  dataset,
  dimred = dynwrap::get_dimred(dataset),
  expression = dataset$expression,
  velocity_vector = dataset$velocity_vector
) {
  assert_that(!is.null(dimred), msg = "Add a dimred to the dataset before adding dimred_future")
  assert_that(
    !is.null(expression),
    !is.null(velocity_vector)
  )

  dataset$dimred_future <- embed_velocity(
    dataset,
    dimred = dimred,
    expression = expression,
    velocity_vector = velocity_vector
  )
  attr(dataset$dimred_future, "dimred_digest") <- digest::digest(dimred, algo = "md5")
  dataset
}


#' @rdname add_dimred_future
#' @export
embed_velocity <- function(
  dataset,
  dimred = dynwrap::get_dimred(dataset),
  expression = dataset$expression,
  velocity_vector = dataset$velocity_vector
) {
  if (check_scvelo(dataset$velocity$scvelo)) {
    dimred_future <- embed_velocity_scvelo(
      dataset$velocity$scvelo,
      dimred = dimred
    )
  } else {
    dimred_future <- embed_velocity_new(
      dataset,
      dimred = dimred,
      expression = expression,
      velocity_vector = velocity_vector
    )
  }

  dimred_future
}


embed_velocity_new <- function(
  dataset,
  dimred = dynwrap::get_dimred(dataset),
  expression = dataset$expression,
  velocity_vector = dataset$velocity_vector
) {
  assert_that(
    !is.null(expression),
    !is.null(velocity_vector),
    !is.null(dimred)
  )

  # create adata object
  velo <- as.matrix(velocity_vector)

  adata = anndata$AnnData(expression)
  adata$var_names <- colnames(expression)
  adata$obs_names <- rownames(expression)

  py_set_item(adata$layers, "spliced", expression)
  py_set_item(adata$layers, "velo", velo)

  # is necessary internally
  adata$uns[["velo_settings"]] <- list(mode = "deterministic")

  # embed velocity
  scvelo$tl$velocity_graph(adata, vkey = "velo")

  # assign dimred and embed
  py_set_item(adata$obsm, "X_dimred", dimred)
  scvelo$tl$velocity_embedding(adata, basis = "dimred", vkey = "velo")

  dimred_future <- dimred + adata$obsm[["velo_dimred"]]
  dimnames(dimred_future) <- dimnames(dimred)

  # replace the ocasional NA dimred position
  dimred_future[is.na(dimred_future)] <- dimred[is.na(dimred_future)]

  dimred_future
}



# embed the velocity directly using an scvelo anndata object, instead of recreating one
embed_velocity_scvelo <- function(adata, dimred) {
  py_set_item(adata$obsm, "X_dimred", dimred)
  scvelo$tl$velocity_embedding(adata, "dimred")
  dimred_future <- dimred + adata$obsm[["velocity_dimred"]]
  dimnames(dimred_future) <- dimnames(dimred)

  # replace the ocasional NA dimred position
  dimred_future[is.na(dimred_future)] <- dimred[is.na(dimred_future)]

  dimred_future
}
