#' Column-centre a matrix  (replicates the centering step in spmDesignOS)
#'
#' @keywords internal

center_cols <- function(M) {
  n       <- nrow(M)
  col_sum <- colSums(M)
  t(apply(M, 1, function(row) row - col_sum / n))
}
