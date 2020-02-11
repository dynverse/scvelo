test_that("velocity can be calculated", {
  set.seed(1)
  unspliced <- Matrix::rsparsematrix(50, 500, 1, rand.x = runif)
  dimnames(unspliced) <- list(paste0('A', seq_len(nrow(unspliced))), paste0('B', seq_len(ncol(unspliced))))
  spliced <- Matrix::rsparsematrix(50, 500, 1, rand.x = runif)
  dimnames(spliced) <- dimnames(unspliced)

  velocity <- get_velocity(spliced, unspliced)

  expect_is(velocity$scvelo, "anndata.core.anndata.AnnData")
})
