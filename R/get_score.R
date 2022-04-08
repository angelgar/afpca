#' Title
#'
#' @param yi
#' @param Mu
#' @param Phi
#' @param sigma
#' @param n.comp
#'
#'
#' @keywords internal
#' @return
#'
#'
#'
#' @examples
get_score <- function(yi, Mu, Phi, sigma, n.comp) {


  score <- Matrix::solve(t(Phi) %*% Phi / sigma + diag(rep(1,n.comp))) %*% t(Phi) %*% (yi - Mu) / sigma
  return(t(score))

}
