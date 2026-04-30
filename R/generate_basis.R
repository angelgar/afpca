#' Generate Spline Basis to estimate FPC
#'
#' @description
#' This is an internal function to generate a spline basis used to calculate the FPCs using asp2 from AdaptFitOS
#'
#' @param data input data
#' @param basis  Type of spline basis. "os" for B-spline basis and "trunc.poly" for truncated polynomials (default)
#' @param poly.degree Degree of spline basis
#' @param knots passed to AdaptFitOS
#' @param nbs number of spline basis
#'
#' @keywords internal
generate_basis <- function(data, basis = "trunc.poly",
                                 poly.degree = 3, knots = NA,
                                 nbs = 40) {


  data <- t(data)


  ## Create dummy data based on the dataset
  X_temp<- (1:dim(data)[2]) / dim(data)[2]
  Y <- rnorm(dim(data)[2])

  if (is.na(knots))
    knots <- default_knots(X_temp, nbs)

  dm         <- spline_design_matrices(X_temp, knots, basis=basis, degree=poly.degree)
  Theta_beta <- dm$Xb
  Theta_b    <- dm$Zb
  Theta      <- dm$Theta

  return(list(Theta = Theta,
              Theta_beta = Theta_beta,
              Theta_b = Theta_b))

}
