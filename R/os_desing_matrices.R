#' O'Sullivan spline decomposition into fixed (C) and random (Z) parts.
#' Mirrors Predict.matrix.lme() / smooth.construct.os.smooth.spec() from asp-internal.r.
#'
#' @param x predictor vector (length n)
#' @param nk number of interior knots (only the COUNT matters; knots are forced equidistant)
#' @param degree length-2 vector c(p, q)  where p = B-spline degree, q = penalty order
#'
#' @keywords internal

os_design_matrices <- function(x, nk, degree) {

  requireNamespace("splines", quietly = TRUE)

  p <- degree[1]
  q <- degree[2]
  n <- length(x)

  # AdaptFitOS always uses equidistant interior knots for 'os', ignoring
  # whatever knot locations were passed – only the count (nk) is used.
  k <- seq(min(x), max(x), length.out = nk + 2)[-c(1, nk + 2)]
  names(k) <- NULL

  # --- B-spline basis (with intercept) ---
  X_bs <- splines::bs(x,
                      knots           = k,
                      degree          = p,
                      Boundary.knots  = c(min(x), max(x)),
                      intercept       = TRUE)
  X_bs <- unclass(X_bs)            # strip 'bs' class so plain matrix arithmetic works
  X_bs <- center_cols(X_bs)       # centre each column

  # --- Build q-th derivative operator d ---
  d <- diag(ncol(X_bs))
  for (i in seq_len(q)) {
    d           <- diff(d)
    aK          <- c(rep(min(x), p + 1 - i), k, rep(max(x), p + 1 - i))
    idx_num     <- seq_len(nrow(d))                      # numerator knot indices
    idx_den_lo  <- seq_len(nrow(d))                      # lower denominator indices
    idx_den_hi  <- idx_den_lo + (p + 1 - i)
    wts         <- 1 / (aK[idx_den_hi] - aK[idx_den_lo])
    d           <- d * (p + 1 - i) * matrix(rep(wts, each = ncol(d)),
                                            ncol = ncol(d),
                                            nrow = nrow(d),
                                            byrow = TRUE)
  }

  # --- Integral of products of derivative B-splines (R_int) ---
  aK_q  <- c(rep(min(x), p - q + 1), k, rep(max(x), p - q + 1))
  sz    <- nk + p + 1 - q
  R_int <- matrix(0, sz, sz)

  for (i in seq_len(sz)) {
    j_max <- min(i + p - q, sz)
    for (j in i:j_max) {
      integrand <- function(u) {
        Nq <- splines::spline.des(aK_q, u,
                                  ord     = p - q + 1,
                                  derivs  = rep(0, length(u)),
                                  outer.ok = TRUE)$design
        Nq[, i] * Nq[, j]
      }
      x1         <- aK_q[j]
      x2         <- aK_q[i + (p - q) + 1]
      R_int[i, j] <- stats::integrate(integrand, x1, x2, subdivisions = 1000L)$value
    }
  }
  R_int         <- R_int + t(R_int)
  diag(R_int)   <- diag(R_int) / 2

  # --- Durban mixed-model decomposition ---
  Re    <- eigen(R_int)
  Re12  <- Re$vectors %*% diag(sqrt(Re$values)) %*% t(Re$vectors)
  D     <- Re12 %*% d
  DI    <- tcrossprod(D)                          # D %*% t(D)

  Z_os  <- X_bs %*% t(D) %*% solve(DI)           # random-effects columns

  # --- Fixed-effects null-space columns (C) ---
  Dq         <- t(d) %*% R_int %*% d
  Oe         <- eigen(Dq)
  null_idx   <- (ncol(X_bs) - q + 1):(ncol(X_bs) - 1)
  U0         <- Oe$vectors[, null_idx, drop = FALSE]
  C_os       <- X_bs %*% U0
  dimnames(C_os)[[2]] <- NULL

  list(C = C_os, Z = Z_os)
}
