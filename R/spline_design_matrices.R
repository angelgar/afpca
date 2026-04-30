#' Main replacement for asp2()'s design matrix extraction.
#' Internal function to generate knots
#'
#' @param x numeric predictor vector (your X_temp)
#' @param knots numeric knot vector (from default_knots() or user-supplied)
#' @param basis "trunc.poly" | "tps" | "os"   (what you pass as 'basis' in f())
#' @param degree integer (or c(p,q) for "os")  (your poly.degree)
#' @keywords internal



spline_design_matrices <- function(x, knots, basis = "os", degree = 3) {

  n <- length(x)

  # ---- "trunc.poly" (truncated power functions) ----------------------------
  if (basis == "trunc.poly") {

    # Fixed effects: intercept + centred x^1, …, x^degree
    Xmat <- matrix(0, n, degree)
    for (d in seq_len(degree)) {
      col          <- x^d
      Xmat[, d]   <- col - sum(col) / n
    }

    # Random effects: centred (x - k_j)^degree_+
    Zmat  <- outer(x, knots, "-")
    Zmat  <- (Zmat * (Zmat > 0))^degree
    Zmat  <- center_cols(Zmat)

    Xb <- cbind(1, Xmat)
    Zb <- Zmat

    # ---- "tps" (thin plate splines) ------------------------------------------
  } else if (basis == "tps") {

    ncol_X <- (degree - 1) / 2
    if (floor(ncol_X) != ncol_X)
      stop("'tps' basis requires an odd degree (e.g., 1, 3, 5).")

    # Fixed effects: intercept + centred x^1, …, x^((degree-1)/2)
    if (ncol_X > 0) {
      Xmat <- matrix(0, n, ncol_X)
      for (d in seq_len(ncol_X)) {
        col        <- x^d
        Xmat[, d] <- col - sum(col) / n
      }
    } else {
      Xmat <- NULL
    }

    # Random effects: TPS basis transformed to canonical mixed-model form
    Zmat        <- abs(outer(x, knots, "-"))^degree
    sqrt_Omega  <- matrix_sqrt(abs(outer(knots, knots, "-"))^degree)
    Zmat        <- t(solve(sqrt_Omega, t(Zmat)))
    Zmat        <- center_cols(Zmat)

    Xb <- cbind(1, Xmat)
    Zb <- Zmat

    # ---- "os" (O'Sullivan penalised B-splines) --------------------------------
  } else if (basis == "os") {

    # Resolve degree to c(p, q) if a scalar was supplied
    if (length(degree) == 1) {
      if (((degree + 1) / 2) %% 1 != 0) {
        warning("'os' basis: degree must satisfy (degree+1)/2 is integer. ",
                "Defaulting to c(3, 2).")
        degree <- c(3L, 2L)
      } else {
        degree <- c(degree, (degree + 1L) / 2L)
      }
    }

    # 'os' uses only the NUMBER of knots (knots are forced equidistant internally)
    nk  <- length(knots)
    cz  <- os_design_matrices(x, nk, degree)

    if (degree[1] != 1L) {
      Xb <- cbind(1, cz$C)
    } else {
      Xb <- matrix(1, n, 1)
    }
    Zb <- cz$Z

  } else {
    stop("Unsupported basis '", basis,
         "'. Choose one of: 'trunc.poly', 'tps', 'os'.")
  }

  colnames(Xb) <- NULL
  colnames(Zb) <- NULL
  Theta <- cbind(Xb, Zb)

  list(Xb = Xb, Zb = Zb, Theta = Theta)
}
