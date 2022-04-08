#' Estimate the Full Likelihood of the FPCA Model (Internal)
#'
#' @param data Input Data i.e. Y
#' @param coef Estimated Spline Coefficients
#' @param Theta Matrix of Penalized Spline Basis
#' @param c_mat Matrix of Scores
#' @param N.Unp.Basis Matrix of Unpenalized Spline Basis
#' @param nbs Number of penalized spline basis
#' @param n.comp Number of FPC's that are being estimated
#' @param lambda Adaptive Smoothness Penalty Matrix
#' @param sigma Residual Variance
#'
#' @return
#'

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
