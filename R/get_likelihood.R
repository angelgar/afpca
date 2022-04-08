#' Title
#'
#' @param data
#' @param coef
#' @param Theta
#' @param c_mat
#' @param N.Unp.Basis
#' @param nbs
#' @param n.comp
#' @param lambda
#' @param sigma
#'
#' @keywords internal
#' @return
#'
#'
#' @examples
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
