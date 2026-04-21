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

  if (is.na(knots)) {
    knots = AdaptFitOS::default.knots(X_temp,nbs)
  }

  ## Run asp2 to get the spline basis

  y.fit <- AdaptFitOS::asp2(Y ~ f(X_temp, basis = basis, degree = poly.degree,
                            adap = F, knots = knots),
                spar.method = "ML")

  Theta_beta = y.fit$design.matrices$Xb
  Theta_b = y.fit$design.matrices$Zb
  Theta <- cbind(Theta_beta, Theta_b)

  return(list(Theta = Theta,
              Theta_beta = Theta_beta,
              Theta_b = Theta_b))

}
