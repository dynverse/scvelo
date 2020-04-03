context("testing get.R")

set.seed(1)
unspliced <- Matrix::rsparsematrix(50, 100, 0.5, rand.x = runif)
dimnames(unspliced) <- list(paste0('A', seq_len(nrow(unspliced))), paste0('B', seq_len(ncol(unspliced))))
spliced <- Matrix::rsparsematrix(50, 100, 0.5, rand.x = runif)
dimnames(spliced) <- dimnames(unspliced)

velocity <- get_velocity(spliced, unspliced)

test_that("get_lambda", {
  lambdas <- get_lambda(velocity)

  expect_true(is.numeric(lambdas))
  expect_true(length(colnames(spliced)) == length(lambdas))
  expect_true(all(colnames(spliced) == names(lambdas)))
})
