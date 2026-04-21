#' Generate Spline Basis to estimate FPC
#'
#' @description
#' This is to generate default knots
#' This code belongs to AdaptFitOS
#' @keywords internal
#'
#'


Predict.matrix.lme <- function (object, data, drv = 0, center = T, ...)
{
  m <- object$m
  q = m[2]
  p = m[1]
  x <- data[[object$term]]
  n = length(x)
  k = object$knots
  nk <- object$bs.dim
  if (inherits(object, "ospline.smooth")) {
    X <- Predict.matrix.ospline.smooth(object, data, drv = drv)
    d = diag(ncol(X))
    for (i in 1:q) {
      d = diff(d)
      allKnots_p <- c(rep(min(x), p + 1 - i), k, rep(max(x),
                                                     p + 1 - i))
      weights = matrix(rep(1/(allKnots_p[-(1:(p + 1 - i))] -
                                allKnots_p[-((nk + 2 * p - p - i + 2):(nk + 2 *
                                                                         p))]), each = ncol(d)), ncol(d), nrow(d))
      d = d * (p + 1 - i) * t(weights)
    }
    allKnots <- c(rep(min(x), (p - q) + 1), k, rep(max(x),
                                                   (p - q) + 1))
    R_int = matrix(0, nk + p + 1 - q, nk + p + 1 - q)
    for (i in 1:(nk + p - q + 1)) {
      for (j in i:(min(i + p - q, nk + p - q + 1))) {
        R <- function(x) {
          Nq <- splines::spline.des(allKnots, x, p - q + 1, derivs = 0 *
                             x, outer.ok = T)$design
          Nq[, i] * Nq[, j]
        }
        x1 = allKnots[j]
        x2 = allKnots[(i + (p - q) + 1)]
        R_int[i, j] <- stats::integrate(R, x1, x2, subdivisions = 1000)$value
      }
    }
    R_int <- R_int + t(R_int)
    diag(R_int) <- diag(R_int)/2
    Re = eigen(R_int)
    Re12 = Re$vectors %*% diag(sqrt(Re$values)) %*% t(Re$vectors)
    D = Re12 %*% d
    DI = tcrossprod(D)
    if (center)
      X = t(apply(X, 1, function(x) x - colSums(X)/n))
    Z = X %*% t(D) %*% solve(DI)
    Dq = t(d) %*% R_int %*% d
    O.e = eigen(Dq)
    null.space = (ncol(X) - q + 1):(ncol(X) - 1)
    U0 = O.e$vectors[, null.space]
    C = X %*% U0
    dimnames(C)[[2]] = NULL
    newknots = seq(min(x), max(x), length = nk + p + 1 -
                     q + 2)[-c(1, nk + p + 1 - q + 2)]
    return(list(C = C, Z = Z, knots = newknots))
  }
  else if (inherits(object, "tlspline.smooth") | inherits(object,
                                                          "trunc.poly")) {
    x <- as.vector(data[[object$term]])
    Z <- outer(x, k, "-")
    Z <- (Z * (Z > 0))^m[1]
    C = rep(1, length(x))
    for (i in 1:(m[1])) C = cbind(C, x^i)
    if (center) {
      n = nrow(C)
      colSC = colSums(C)
      colSZ = colSums(Z)
      C = t(apply(C, 1, function(x) x - colSC/n))
      Z = t(apply(Z, 1, function(x) x - colSZ/n))
    }
  }
  else if (inherits(object, "tps") | inherits(object, "ts.smooth")) {
    if (is.null(m))
      m = object$p.order
    if (is.null(k))
      stop("No knots given in smooth.construct.")
    x <- as.vector(data[[object$term]])
    svd.Omega = svd(abs(outer(k, k, "-"))^m[1])
    matrix.sqrt.Omega = t(svd.Omega$v %*% (t(svd.Omega$u) *
                                             sqrt(svd.Omega$d)))
    Z = t(solve(matrix.sqrt.Omega, t(abs(outer(x, k, "-")^m[1]))))
    C = cbind(rep(1, length(x)))
    for (i in 1:((m[1] - 1)/2)) C = cbind(C, x^i)
    if (center) {
      n = nrow(C)
      colSC = colSums(C)
      colSZ = colSums(Z)
      C = t(apply(C, 1, function(x) x - colSC/n))
      Z = t(apply(Z, 1, function(x) x - colSZ/n))
    }
  }
  else stop("scbM can be fitted only with os, tl or tps basis functions")
  list(C = C, Z = Z)
}
