#' Direct replacement for AdaptFitOS::default.knots() / SemiPar::default.knots()
#' Internal function to generate knots
#'
#' @param x numeric vector (the predictor)
#' @param num.knots integer, number of knots; auto-computed when omitted
#' @param knotchoice "quantiles" (default) or "equidistant"
#' @keywords internal


default_knots <- function(x, num.knots, knotchoice = "quantiles") {

  x <- unique(x)

  if (missing(num.knots)) {
    n         <- length(x)
    d         <- max(4, floor(n / 35))
    num.knots <- floor(n / d - 1)
  }

  x <- as.vector(x[!is.na(x)])

  if (knotchoice == "equidistant") {
    knots <- seq(min(x), max(x), length.out = num.knots)
  } else if (knotchoice == "quantiles") {
    probs <- seq(0, 1, length.out = num.knots + 2)[-c(1, num.knots + 2)]
    knots <- stats::quantile(x, probs, names = FALSE)
  } else {
    stop("Unknown knotchoice. Use 'quantiles' or 'equidistant'.")
  }

  names(knots) <- NULL
  knots
}
