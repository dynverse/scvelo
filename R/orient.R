#' Reorients the edges of the milestone network to the cell's RNA velocity vectors
#'
#' @inheritParams dynwrap::common_param
#' @inheritParams dynwrap::add_expression
#'
#' @return The trajectory with oriented *milestone_network* and *progressions*
#'
#' @examples
#' # we'll create a simple linear trajectory
#' cell_ids <- c("a", "b", "c", "d", "e")
#' pseudotime <- setNames(seq_along(cell_ids), cell_ids)
#' expression <- as.matrix(data.frame(
#'   a = pseudotime,
#'   b = pseudotime ** 2,
#'   c = log(pseudotime)
#' ))
#' expression_future <- as.matrix(data.frame(
#'   a = (pseudotime + 1),
#'   b = (pseudotime + 1) ** 2,
#'   c = log(pseudotime + 1)
#' ))
#'
#' # the milestone network is "wrong" in the sense that B and A are oriented in the opposite direction
#' milestone_network <- tibble::tribble(
#'   ~from, ~to, ~length, ~directed,
#'   "B", "A", 1, TRUE,
#'   "B", "C", 1, TRUE
#' )
#' progressions <- tibble::tribble(
#'   ~cell_id, ~from, ~to, ~percentage,
#'   "a", "B", "A", 1,
#'   "b", "B", "A", 0.5,
#'   "c", "B", "A", 0,
#'   "d", "B", "C", 0.5,
#'   "e", "B", "C", 1
#' )
#'
#' trajectory <- dynwrap::wrap_expression(
#'   counts = expression,
#'   expression = expression,
#'   expression_future = expression_future
#' )
#'
#' trajectory <- dynwrap::add_trajectory(
#'   trajectory,
#'   milestone_network = milestone_network,
#'   progressions = progressions
#' )
#'
#' trajectory_oriented <- orient_topology_to_velocity(trajectory)
#'
#' # the edge is now correctly oriented
#' trajectory_oriented$milestone_network
#' assertthat::assert_that(
#'   all(
#'     trajectory_oriented$milestone_network[2, c("from", "to")] == c("B", "C")
#'   )
#' )
#'
#' @export
orient_topology_to_velocity <- function(
  trajectory,
  velocity = trajectory$velocity,
  expression_source = trajectory,
  expression = dynwrap::get_expression(expression_source, "expression"),
  expression_future = dynwrap::get_expression(expression_source, "expression_future")
) {
  # dummy proofing
  assert_that(is(trajectory, "dynwrap::with_trajectory"))
  assert_that(!is.null(expression))
  assert_that(!is.null(expression_future))

  if (nrow(trajectory$divergence_regions)) {
    stop("Orienting topologies with divergence regions doesn't work yet")
  }

  # flip_fractions <- pmap_dbl(trajectory$milestone_network, function(from, to, ...) {
  #   # order cells based on their percentage
  #   progressions_edge <- trajectory$progressions %>%
  #     filter(from == !!from, to == !!to) %>%
  #     arrange(desc(percentage))
  #
  #   # find for each cell its nearest neighbor (not self) in the expression_future
  #   nn_ix <- FNN::knnx.index(
  #     expression[progressions_edge$cell_id, ],
  #     expression_future[progressions_edge$cell_id, ],
  #     k = 2
  #   )
  #   nn_ix[, 1][nn_ix[, 1] == seq_len(nrow(nn_ix))] <- NA
  #   nn_ix <- nn_ix %>% apply(1, function(x) first(x[!is.na(x)]))
  #
  #   # find out whether RNA velocity support and edge reversal
  #   sum(nn_ix > seq_along(nn_ix))/sum(nn_ix < seq_along(nn_ix))
  # })
  #
  # milestone_network_toflip <- trajectory$milestone_network %>%
  #   filter(flip_fractions > 1)
  #
  # trajectory <- flip_edges(trajectory, milestone_network_toflip)

  assert_that(!is.null(velocity$transition_matrix))
  transition_matrix <- velocity$transition_matrix
  assert_that(nrow(transition_matrix) == nrow(expression))
  assert_that(ncol(transition_matrix) == nrow(expression))

  flip_fractions <- pmap_dbl(trajectory$milestone_network, function(from, to, ...) {
    progressions_edge <- trajectory$progressions %>%
      filter(from == !!from, to == !!to) %>%
      arrange(percentage)

    transition_matrix_ordered <- Matrix::t(transition_matrix[progressions_edge$cell_id, progressions_edge$cell_id])

    sum(Matrix::tril(transition_matrix_ordered, 1)) / sum(Matrix::triu(transition_matrix_ordered, 1))
  })
  trajectory$milestone_network %>% mutate(fraction = flip_fractions)
  trajectory <- flip_edges(trajectory, trajectory$milestone_network %>% mutate(fraction = flip_fractions) %>% filter(fraction < 1))

  # remove the root
  trajectory <- trajectory %>% remove_root()

  # tada!
  trajectory
}

