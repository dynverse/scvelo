get_lambda <- function(velocity) {
  lambdas <- velocity$scvelo$var[["velocity_gamma"]]
  names(lambdas) <- velocity$scvelo$var_names$to_list()
  lambdas
}
