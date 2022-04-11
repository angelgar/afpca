#' Simulate Adaptive Functional Data
#'
#'
#' @description
#' Simulated functional data that includes with sharp changes over the domain of the function (i.e. smoothness rapidly changes over the domain). The data simulates data using sinusoidal functions with varying period and amplitude.
#'
#' @param n.tp number of timepoints
#' @param N.subj number of subjects
#' @param var.scores variance of the scores (this is the value of the eigenvalues)
#' @param noise.var residual noise variance
#' @param seed.num seed number
#'
#' @importFrom magrittr %>%
#'
#'
#' @return An object of class "fpca_sim_data" that contains the following:
#'
#' \item{data}{simulated dataset}
#' \item{data_true}{simulated dataset without any noise}
#' \item{Mu_true}{The true data-generating mean-function}
#' \item{Phi_true}{The true data-generating functional principal components}
#' \item{Scores_true}{The true scores}
#'
#' @export
#'
#' @examples
#' sim_data <- simulate_adaptive_functional_data(n.tp = 100, N.subj = 25)
simulate_adaptive_functional_data <- function(n.tp = 200,
                                              N.subj = 20,
                                              var.scores = c(4,1),
                                              noise.var = 0.05,
                                              seed.num = 1) {


  set.seed(seed.num)

  n.tp.sin <- n.tp / 2
  t = seq(.0625, 1, length.out = n.tp.sin)^(1/4)
  tprime = t^(-3/2)
  period <- 4 * pi * t

  Mu_true <- c(rep(0,n.tp.sin), -rev(tprime * sin((1/2)*period))) * 0.2

  ## Sine function with varying period
  Phi_true <- cbind(c(rep(0,n.tp.sin), -rev(tprime * sin(2*period))),
                    c(rep(0,n.tp.sin), -rev(tprime * sin(1*period))))

  ## Check that it's orthogonormal
  t(Phi_true) %*% Phi_true

  ## I'm just forcing it so that it is numerically orthonormal
  Phi_true <- Phi_true %>% matlib::GramSchmidt(normalize = T)

  Phi_true <- Phi_true*-1

  ## True Scores

  Scores_true <- MASS::mvrnorm(n = N.subj, rep(0,2),
                               Sigma = diag(var.scores),
                               empirical = T)

  ## Generate Data

  data_true <- matrix(NA, nrow = n.tp, ncol = N.subj)

  for (i in 1:N.subj) {

      data_true[,i] <- Mu_true + Phi_true %*% Scores_true[i, ]

  }

  data <- matrix(NA, nrow = n.tp, ncol = N.subj)


  for (i in 1:N.subj) {

    data[,i] <- Mu_true + Phi_true %*% Scores_true[i, ] + rnorm(n.tp, mean = 0, sd = noise.var)

  }

  output <- list(data = data,
                 data_true = data_true,
                 Mu_true = Mu_true,
                 Phi_true = Phi_true,
                 Scores_true = Scores_true)

  class(output) <- "fpca_sim_data"

  return(output)

}
