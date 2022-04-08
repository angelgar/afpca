get_score <- function(yi, Mu, Phi, sigma, n.comp) {


  score <- Matrix::solve(t(Phi) %*% Phi / sigma + diag(rep(1,n.comp))) %*% t(Phi) %*% (yi - Mu) / sigma
  return(t(score))

}
