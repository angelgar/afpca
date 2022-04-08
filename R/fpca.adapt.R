#' Adpative Functional Principal Component Analysis
#'
#' @param data
#' @param poly.degree
#' @param knots
#' @param nbs
#' @param n.comp
#' @param normalize.scores
#' @param orthogonalize_scores
#' @param seed.num
#' @param orthogonalize_fpcs
#' @param ntimes
#' @param basis
#' @param pve
#'
#' @return
#' @export
#'
#' @examples
#'
#'
fpca.adapt <- function(data, poly.degree = 2,
                           knots = NA,
                           nbs = 35,
                           n.comp = 18,
                           normalize.scores = TRUE,
                           orthogonalize_scores = TRUE,
                           seed.num = 4,
                           orthogonalize_fpcs = TRUE,
                           ntimes = 100,
                           basis = "trunc.poly", pve = 0.99) {

  set.seed(seed.num)

  ## Generate Spline Basis
  Basis <- generate_basis(data, poly.degree, knots,
                             nbs,
                             basis)
  Theta <- Basis$Theta
  N.Unp.Basis <- dim(Basis$Theta_beta)[2]
  nbs <- dim(Basis$Theta_b)[2]

  ## Create Initial Values (Random Scores / OLS Estimates)
  Init_Estimates <- init_par(data, n.comp,
                             Theta,
                             orthogonalize_fpcs,
                             N.Unp.Basis)

  N.subj <- dim(data)[2]

  ## Estimate Scores Given Initial Values
  c_mat <- foreach (i=1:N.subj, .combine = rbind) %do% {
    get_score(data[,i], Mu = Init_Estimates$Mu_init, Phi = Init_Estimates$Phi_init,
              sigma = Init_Estimates$sigma_init, n.comp = n.comp)
  }

  c_mat <- orthonormal_step(c_mat, orthogonalize_scores,
                            normalize.scores)


  ## Get Initial Adaptive FPCA Estimates
  coef <- get_coefficients(data = data, Theta,
                           c_mat = c_mat, Lambda = Init_Estimates$lambda_init,
                           sigma = Init_Estimates$sigma_init,
                           N.Unp.Basis, nbs)


  est_fpc <- get_fpcs(Theta, coef,
                      nbs, n.comp, orthogonalize_fpcs,
                      N.Unp.Basis)


  ## Get Sigma
  sigma <- get_sigma(data, coef, Theta, c_mat,
                     N.Unp.Basis, nbs)

  lambda <- get_lambda(coef, n.comp,
                       N.Unp.Basis)

  loglik <- get_likelihood(data, coef, Theta, c_mat,
                           N.Unp.Basis, nbs, n.comp, lambda, sigma)


  ## Iterate until convergence

  no.iter <- 1

  for (times in 1:ntimes) {

    ## Update Scores
    c_mat <- foreach (i=1:N.subj,.combine = rbind) %do% {
      get_score(data[,i], Mu = est_fpc$Mu, Phi = est_fpc$Phi,
                sigma = sigma, n.comp = n.comp)
    }

    c_mat <- orthonormal_step(c_mat, orthogonalize_scores,
                              normalize.scores)

    ## Update Lambda
    lambda <- get_lambda(coef, n.comp,
                         N.Unp.Basis)


    # Update FPCs

    coef_old <- coef

    coef <- get_coefficients(data = data, Theta,
                             c_mat = c_mat, Lambda = lambda,
                             sigma = sigma, N.Unp.Basis, nbs)


    est_fpc <- get_fpcs(Theta, coef,
                        nbs, n.comp, orthogonalize_fpcs, N.Unp.Basis)

    sigma_old <- sigma
    sigma <- get_sigma(data, coef, Theta, c_mat, N.Unp.Basis, nbs)

    loglik_old <- loglik
    loglik <- get_likelihood(data, coef, Theta, c_mat,
                             N.Unp.Basis, nbs, n.comp, lambda, sigma)



    ## Check for Convergence
    # if (cor(coef, coef_old) > 0.999 & times > 4) {
    #  break
    #}

    no.iter <- no.iter + 1
    print(loglik)
    print(loglik_old)
    (loglik - loglik_old)/ abs(loglik_old)


    print(no.iter)
    print(sigma_old)
    print(sigma)
    convergence <- FALSE
    if (((loglik - loglik_old)/ abs(loglik_old)) < 0.0001 & times > 4) {
      convergence <- TRUE
      break
    }
  }

  est_fpc <- get_fpcs(Theta, coef,
                      nbs, n.comp, orthogonalize_fpcs = T, N.Unp.Basis)

  eigen.values <- est_fpc$eigen.vals
  fpc_no_keep <- min(which(cumsum(eigen.values^2) / sum(eigen.values^2) > pve))

  ## Update Scores
  c_mat <- foreach (i=1:N.subj,.combine = rbind) %do% {
    get_score(data[,i], Mu = est_fpc$Mu, Phi = est_fpc$Phi[, 1:fpc_no_keep],
              sigma = sigma, n.comp = fpc_no_keep)
  }

  reconstructed_data <- matrix(NA, nrow = dim(data)[1], ncol = dim(data)[2])

  for (i in 1:dim(data)[2]) {

    # reconstructed_data[,i] <- est_fpc$Mu + est_fpc$Phi[,] %*% c_mat[i,]
    reconstructed_data[,i] <- est_fpc$Mu + est_fpc$Phi[,1:fpc_no_keep] %*% c_mat[i,]

  }

  return(list(Y = data,
              Y_hat = reconstructed_data,
              mean = est_fpc$Mu,
              fpcs = est_fpc$Phi[,1:fpc_no_keep],
              scores = c_mat,
              eigen.values = est_fpc$eigen.vals,
              Theta = Theta,
              Basis = Basis,
              coef = coef,
              coef_old = coef_old,
              convergence = convergence,
              no.iter = no.iter))

}
