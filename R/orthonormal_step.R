#' Function to orthonormalize (or orthogonalize) scores
#'
#' @description
#' Internal function to orthonormalize or orthogonalize scores
#'
#' @param c_mat matrix of scores
#' @param orthogonalize_scores logical; wether to orthogonalize scores
#' @param normalize.scores logical; wether to normalize scores
#'
#'
#' @keywords internal
orthonormal_step <- function(c_mat, orthogonalize_scores,
                             normalize.scores) {


  if (orthogonalize_scores) {
    if (normalize.scores) {
      c_mat <- qr.Q(qr(c_mat)) %>%
        scale()
    } else {
      c_mat <- qr.Q(qr(c_mat))
    }
  }

  return(c_mat)

}
