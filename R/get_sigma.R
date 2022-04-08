get_sigma <- function(data, coef, Theta, c_mat, N.Unp.Basis, nbs) {


  Y <- matrix(data, ncol = 1)
  N.subj <- dim(data)[2]
  n.comp <- dim(c_mat)[2]

  ## Make The Design Matrix
  C_Big <- cbind(rep(1, N.subj), c_mat)
  Complete_Basis_noshift <- C_Big %x% Theta
  Complete_Basis <- Complete_Basis_noshift %>% as_tibble() %>%
    relocate(as.vector(outer(1:N.Unp.Basis, 1:n.comp * (nbs + N.Unp.Basis), FUN = "+")), .after = 1:N.Unp.Basis) %>%
    as.matrix()

  sigma <- sum((Y - Complete_Basis %*% coef)^2) / length(Y)
  return(sigma)

}
