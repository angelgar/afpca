#' Estimate the Adaptive Penalty Matrix (Lambda)
#'
#' @param coef Coefficients
#' @param n.comp Number of FPCs to estimate
#' @param N.Unp.Basis Number of Unpenalized Basis
#'
#' @keywords internal
#' @return
#'
#' @examples
#'
get_lambda <- function(coef, n.comp, N.Unp.Basis) {


  beta_vec <- coef[1:(N.Unp.Basis*(1 + n.comp))]
  b_vec <- coef[-(1:(N.Unp.Basis*(1 + n.comp)))]

  unpenalized_lambda <- matrix(0, nrow = N.Unp.Basis*(1 + n.comp), ncol = N.Unp.Basis*(1 + n.comp))

  penalized_lambda <- diag(ifelse(b_vec^(-2) > 1e+7, 1e+7, b_vec^(-2)))

  lambda <- Matrix::bdiag(unpenalized_lambda, penalized_lambda) %>%
    as.matrix()

  return(lambda)

}
