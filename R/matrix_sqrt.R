#' Positive-definite matrix square root via SVD  (from SemiPar::matrix.sqrt)
#'
#' @keywords internal


matrix_sqrt <- function(A) {
  sva <- svd(A)
  if (min(sva$d) < 0) stop("Matrix square root is not defined (negative eigenvalue).")
  t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
}
