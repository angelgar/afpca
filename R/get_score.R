#' Generate subject-specific scores given estimates of mean and eigenfunctions
#'
#' @description
#' Internal function to calculate scores
#'
#' @param yi one function input data
#' @param Mu estimated mean
#' @param Phi matrix of estimated fpcs
#' @param sigma residual variance
#' @param n.comp number of estimated fpcs
#'
#'
#' @keywords internal
#' @return
#'
get_score <- function(yi, Mu, Phi, sigma, n.comp) {


  score <- Matrix::solve(t(Phi) %*% Phi / sigma + diag(rep(1,n.comp))) %*% t(Phi) %*% (yi - Mu) / sigma
  return(t(score))

}
