# Internal utility functions for the fql package
# These are not exported

#' Generalized Inverse of a Matrix
#'
#' Computes the Moore-Penrose generalized inverse using SVD.
#'
#' @param X A numeric matrix.
#' @param tol Tolerance for singular value thresholding.
#' @return The generalized inverse matrix.
#' @keywords internal
.ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  dnx <- dimnames(X)
  if (is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if (any(nz)) s$v[, nz] %*% (t(s$u[, nz]) / s$d[nz]) else X,
    dimnames = dnx[2:1]
  )
}

#' Kernel Variance Estimation
#'
#' Estimates the variance at a single point using a Gaussian kernel smoother
#' (local linear regression of squared residuals on fitted means).
#'
#' @param mi The point at which to estimate variance.
#' @param mu Vector of all fitted means.
#' @param er Vector of squared residuals (y - mu)^2.
#' @param bw Bandwidth for the Gaussian kernel.
#' @return Estimated variance at point mi.
#' @keywords internal
.kernel_varest <- function(mi, mu, er, bw) {
  n <- length(mu)
  ker <- exp(-((mu - mi) / bw)^2 / 2)
  An0 <- sum(ker) / (n * bw)
  An1 <- sum(ker * (mu - mi)) / (n * bw)
  An2 <- sum(ker * (mu - mi)^2) / (n * bw)
  WI <- (ker * (An2 - (mu - mi) * An1) / (An0 * An2 - An1^2)) / (n * bw)
  sum(WI * er)
}

#' Compute Linear Predictor
#'
#' Computes eta = X %*% beta, optionally adding an offset.
#'
#' @param x Design matrix (n x p).
#' @param beta Coefficient vector (length p).
#' @param offset Optional offset vector (length n) or NULL.
#' @return Linear predictor vector eta (length n).
#' @keywords internal
.compute_eta <- function(x, beta, offset = NULL) {
  eta <- as.vector(x %*% beta)
  if (!is.null(offset)) {
    eta <- eta + offset
  }
  eta
}

#' Apply Link Function (inverse)
#'
#' Applies the inverse link function to convert from linear predictor
#' to mean scale. Currently supports log link only.
#'
#' @param eta Linear predictor vector.
#' @return Mean vector mu.
#' @keywords internal
.link_inv <- function(eta) {
  exp(eta)
}

#' Compute D_i Matrix (Vectorized)
#'
#' For log link: D_i = mu_i * x_i, so Di (n x p) = diag(mu) %*% X.
#' This computes the full Di matrix using vectorized operations.
#'
#' @param mu Fitted mean vector (length n).
#' @param x Design matrix (n x p).
#' @return Di matrix (n x p) where row i is mu[i] * x[i,].
#' @keywords internal
.compute_Di <- function(mu, x) {
  mu * x  # recycling: each row of x multiplied by corresponding mu
}

#' Safe Matrix Solve
#'
#' Solves a linear system, falling back to generalized inverse if
#' the matrix is rank-deficient.
#'
#' @param A Square matrix.
#' @param b Right-hand side vector.
#' @param p Expected full rank.
#' @return Solution vector.
#' @keywords internal
.safe_solve <- function(A, b, p) {
  r <- qr(A)$rank
  if (r == p) {
    qr.solve(A, b)
  } else {
    .ginv(A) %*% b
  }
}

#' Check Matrix Rank
#'
#' @param A A matrix.
#' @return The rank of A.
#' @keywords internal
.matrix_rank <- function(A) {
  qr(A)$rank
}
