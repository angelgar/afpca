get_fpcs <- function(Theta,
                     coef, nbs, n.comp,
                     orthogonalize_fpcs,
                     N.Unp.Basis) {



  unpenalized.index <- as.vector(outer(1:N.Unp.Basis, 0:n.comp * (nbs + N.Unp.Basis), FUN = "+"))
  Mu <- Theta %*% coef[c(1:N.Unp.Basis, 1:nbs + length(unpenalized.index))]
  Phi <- matrix(NA, nrow = dim(Theta)[1], ncol = n.comp)

  for (index in 1:n.comp) {
    Phi[,index] <- Theta %*% coef[c(1:N.Unp.Basis + N.Unp.Basis*index, 1:nbs + length(unpenalized.index) + nbs*index)]
  }

  if (orthogonalize_fpcs) {
    eigen.vals <- svd(Phi)$d
    Phi <- svd(Phi)$u
  } else {
    eigen.vals <- rep(1,n.comp)
  }

  return(list(Mu = Mu,
              Phi = Phi,
              eigen.vals = eigen.vals))

}
