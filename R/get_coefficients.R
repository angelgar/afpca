#' Get the spline basis coefficients
#'
#' @description
#' Internal function to estimate the spline basis coefficients given the scores and adaptive penalty
#'
#' @param data data input data
#' @param c_mat matrix of scores
#' @param Theta spline basis
#' @param Lambda penalty matrix
#' @param sigma residual variance
#' @param N.Unp.Basis number of unpelized spline basis
#' @param nbs total number of spline basis
#'
#'
#' @keywords internal
get_coefficients <- function(data, c_mat, Theta, Lambda, sigma,
                             N.Unp.Basis, nbs) {


  t <- dim(data)[1]
  N.subj <- dim(data)[2]
  Y <- matrix(data, ncol = 1)
  n.comp <- dim(c_mat)[2]

  ## Create Big Spline Basis and Reorder Coefficients to match the penalty matrix (no penalty for some)
  C_Big <- cbind(rep(1, N.subj), c_mat)
  Complete_Basis_noshift <- C_Big %x% Theta
  Complete_Basis <- Complete_Basis_noshift %>%
                    tibble::as_tibble(.name_repair = ~ vctrs::vec_as_names(..., repair = "unique", quiet = TRUE)) %>%
                    dplyr::relocate(as.vector(outer(1:N.Unp.Basis, 1:n.comp * (nbs + N.Unp.Basis), FUN = "+")), .after = 1:N.Unp.Basis) %>%
                    as.matrix()

  ## Get Coefficients
  coeff <- Matrix::solve(t(Complete_Basis) %*% Complete_Basis / sigma + Lambda, tol = 1e-20) %*% t(Complete_Basis) %*% (Y) / sigma
  return(coeff)

}
