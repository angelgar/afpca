#' Generate Spline Basis to estimate FPC
#'
#' @description
#' This is to generate default knots
#' This code belongs to AdaptFitOS
#' @keywords internal
#'
#'

Predict.matrix.ospline.smooth <- function (object, data, drv = 0, ...)
{
  X <- bs2(data[[object$term]], knots = object$knots, degree = object$m[1],
           Boundary.knots = c(min(data[[object$term]]), max(data[[object$term]])),
           intercept = T, drv = drv)
  X
}
