init_par <- function(data, n.comp,
                     Theta,
                     orthogonalize_fpcs,
                     N.Unp.Basis) {


  t <- dim(data)[1]
  N.subj <- dim(data)[2]
  Y <- matrix(data, ncol = 1)

  ## I'm not sure anymore but should this not be empirical?
  ## I think this is a simple way to enforce the scores being uncorrelated at the first pass

  # scores_init <- mvrnorm(n = N.subj, rep(0,n.comp),
  #                       Sigma = diag(rep(1,n.comp)),
  #                       empirical = T) %>% scale() %>%
  #  as.matrix()

  scores_init <- fpca.face(t(data))$scores %>%
    matlib::GramSchmidt(normalize = FALSE) %>%
    scale() %>%
    as.matrix()

  n.scores.face <- dim(scores_init)[2]

  if (n.scores.face < n.comp) {

    scores_init <- cbind(scores_init,
                         mvrnorm(n = N.subj, rep(0,n.comp - n.scores.face),
                                 Sigma = diag(rep(1,n.comp - n.scores.face)),
                                 empirical = T)) %>%
      matlib::GramSchmidt(normalize = FALSE) %>%
      scale() %>%
      as.matrix()

  } else {

    scores_init <- scores_init[,1:n.comp]

  }



  ## Define the C_init Matrix Vector of 1 and Scores to use in Kronecker

  C_init <- cbind(rep(1, N.subj), scores_init)


  ## Define the Big Spline Basis
  Theta_init <- C_init %x% Theta

  ## Get the OLS Coefficients of the Spline Basis
  lm_init <- lm(Y ~ Theta_init - 1)
  coef_init <- coef(lm_init)

  ## Generate Initial values for Mu and fpcs

  coef.length <- dim(Theta)[2]
  Mu_init <- Theta %*% coef_init[c(1:coef.length)]

  Phi_init <- matrix(NA, nrow = t, ncol = n.comp)
  for (index in 1:n.comp) {
    Phi_init[,index] <- Theta %*% coef_init[c(1:coef.length) + coef.length*index]
  }

  ## Orthonormalize fpcs if needed
  if (orthogonalize_fpcs) {
    Phi_init <- svd(Phi_init)$u
  }

  ## Generate Initial Sigma
  sigma_init <- sum((Y - Theta_init %*% coef_init)^2) / length(Y)

  ## Generate Initial Penalty Matrix
  lambda_init <- get_lambda(coef = coef_init, n.comp = n.comp,
                            N.Unp.Basis)

  return(list(Mu_init = Mu_init,
              Phi_init = Phi_init,
              sigma_init = sigma_init,
              lambda_init = lambda_init))

}
