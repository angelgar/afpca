#' Title
#'
#' @param c_mat
#' @param orthogonalize_scores
#' @param normalize.scores
#'
#'
#' @keywords internal
#' @return
#'
#'
#'
#' @examples
orthonormal_step <- function(c_mat, orthogonalize_scores,
                             normalize.scores) {


  if (orthogonalize_scores) {
    if (normalize.scores) {
      c_mat <- c_mat %>%
        matlib::GramSchmidt(normalize = FALSE) %>%
        scale()
    } else {
      c_mat <- c_mat %>%
        matlib::GramSchmidt(normalize = FALSE)
    }
  }

  return(c_mat)

}
