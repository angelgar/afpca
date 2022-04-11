#' Title
#'
#' @param data
#' @param poly.degree
#' @param knots
#' @param nbs
#' @param basis
#'
#' @keywords internal
#' @return
#'
#' @examples
generate_basis <- function(data, poly.degree = 2, knots = NA,
                              nbs = 40,
                              basis = "trunc.poly") {


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
