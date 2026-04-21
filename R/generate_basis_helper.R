#' Generate Spline Basis to estimate FPC
#'
#' @description
#' This is to generate spline basis
#' This code belongs to AdaptFitOS
#' @keywords internal

generate_basis_helper <- function (form, spar.method = "REML", contrasts = NULL, omit.missing = NULL,
                                   returnFit = FALSE, niter = 20, niter.var = 50, tol = 1e-06,
                                   tol.theta = 1e-06, control = NULL)
  {
  epsilon.fit = epsilon.theta = NULL
  spm.info <- aspmFormReadOS(form, omit.missing, constrasts = contrasts)
  design.info <- spmDesignOS(spm.info)
  Xb <- design.info$X
  Zb <- design.info$Z
  Wb <- cbind(Xb, Zb)
  asp.info <- NULL
  model.matrices <- list(Xb = Xb, Zb = Zb, Wb = Wb)

  asp.info <- list(design.matrices = model.matrices)
  asp.fit.obj <- asp.info
  return(asp.fit.obj)

}
