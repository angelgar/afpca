#' Estimate likelihood of the data
#'
#' @description
#' Internal function to estimate the likelihood of the data, used to check algorithm converges
#'
#' @param data input data
#' @param coef spline coefficients
#' @param Theta spline basis
#' @param c_mat matrix of scores
#' @param N.Unp.Basis number of unpenalized spline basis
#' @param nbs total number of spline basis
#' @param n.comp number of components being estimated
#' @param lambda penalty function
#' @param sigma residual variance
#'
#' @keywords internal
get_likelihood <- function(data, coef, Theta, c_mat,
                           N.Unp.Basis, nbs, n.comp, lambda, sigma) {

  beta_vec <- coef[1:(N.Unp.Basis*(1 + n.comp))]
  b_vec <- coef[-(1:(N.Unp.Basis*(1 + n.comp)))]

  penalty <- diag(lambda)[diag(lambda) > 0 ]

  beta_term <- t(b_vec) %*% diag(penalty) %*% b_vec
  sigma_term <- 0.5*log(sigma)
  lambda_det_term <- 0.5 * sum(log((penalty)))
  scores_term <- sum(c_mat^2)

  full_likelihood <- sigma_term + scores_term + lambda_det_term + beta_term


}
