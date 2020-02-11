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
  assertthat::assert_that(!is.null(velocity$expression_future))

  dataset$expression_future <- velocity$expression_future

  velocity$expression_future <- NULL
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
get_velocity <- function(spliced, unspliced, mode = "stochastic", n_neighbors = 20L) {
  assertthat::assert_that(all(dim(spliced) == dim(unspliced)))
  assertthat::assert_that(all(rownames(spliced) == rownames(unspliced)))
  assertthat::assert_that(all(colnames(spliced) == colnames(unspliced)))

  # filter features with not enough counts in either matrix
  # otherwise scvelo will error quite often (e.g. when calculating the dynamics)
  # feature_ixs <- which((apply(spliced, 2, sd) > 0) & (apply(unspliced, 2, sd) > 0) & (apply(spliced - unspliced, 2, sum) != 0))
  feature_ixs <- rep(TRUE, ncol(spliced))
  spliced <- spliced[, feature_ixs]
  unspliced <- unspliced[, feature_ixs]

  # create anndata object
  velocity = anndata$AnnData(spliced)
  velocity$var_names <- colnames(spliced)
  velocity$obs_names <- rownames(spliced)

  py_assign(velocity$layers, "spliced", spliced)
  py_assign(velocity$layers, "unspliced", unspliced)

  # calculate velocity
  # py_capture_output({ # can't capture output because of https://github.com/rstudio/reticulate/issues/386, otherwise crash when testing

    scvelo$pp$moments(velocity, n_neighbors = n_neighbors)
    if (mode == "dynamical") {
      scvelo$tl$velocity(velocity, mode="deterministic")
      scvelo$tl$velocity_graph(velocity)

      scvelo$tl$recover_dynamics(velocity)
    }
    scvelo$tl$velocity(velocity, mode=mode)
    scvelo$tl$velocity_graph(velocity)
  # })

  velocity_vector <- velocity$layers[["velocity"]]
  velocity_vector[is.na(velocity_vector)] <- 0
  expression_future <- spliced + velocity_vector

  expression_future <- as(expression_future, "dgCMatrix")

  # get transition matrix
  py$x <- velocity
  py_run_string("import scvelo")
  transition_matrix <- py_to_r(py_eval("scvelo.tl.transition_matrix(x).tocsc()"))
  colnames(transition_matrix) <- rownames(spliced)
  rownames(transition_matrix) <- rownames(spliced)

  # return velocity object
  tibble::lst(
    expression_future = expression_future,
    transition_matrix = transition_matrix,
    scvelo = velocity
  )
}


check_scvelo <- function(scvelo) {
  if(is.null(scvelo)) {
    FALSE
  } else if(as.character(scvelo) == "<pointer: 0x0>") {
    FALSE
  } else {
    TRUE
  }
}
