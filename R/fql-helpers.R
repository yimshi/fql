# Core algorithm helper functions for FQL estimation
# These implement the iterative estimation procedure from Shi et al. (2023)

#' IRLS Update for Beta Estimation
#'
#' Performs a single Iteratively Reweighted Least Squares (IRLS) step to
#' update beta, given current variance estimates. Implements Equation (6)
#' from the paper using vectorized matrix operations.
#'
#' @param x Design matrix (n x p).
#' @param y Response vector (length n).
#' @param beta_old Current beta estimate (length p).
#' @param v Variance vector (length n). If NULL, assumes constant variance.
#' @param offset Optional offset vector.
#' @return Updated beta vector (length p).
#' @keywords internal
.irls_step <- function(x, y, beta_old, v = NULL, offset = NULL) {
  p <- ncol(x)
  eta <- .compute_eta(x, beta_old, offset)
  mu <- .link_inv(eta)
  Di <- .compute_Di(mu, x)  # n x p

  if (is.null(v)) {
    # Constant variance (v_i = 1): DV = sum Di Di^T
    # Vectorized: DV = t(Di) %*% Di = crossprod(Di)
    DV <- crossprod(Di)
    DY <- crossprod(Di, y - mu)
  } else {
    # Weighted by 1/v_i: DV = sum Di Di^T / v_i
    # Vectorized: DV = crossprod(Di / sqrt(v))
    w <- 1 / sqrt(as.vector(v))
    Di_w <- Di * w
    DV <- crossprod(Di_w)
    DY <- crossprod(Di, (y - mu) / as.vector(v))
  }

  beta_new <- beta_old + as.vector(.safe_solve(DV, DY, p))
  beta_new
}

#' IRLS Loop for Beta Estimation
#'
#' Runs the full IRLS iteration until convergence to estimate beta,
#' given fixed variance estimates.
#'
#' @param x Design matrix (n x p).
#' @param y Response vector (length n).
#' @param beta_init Initial beta estimate (length p).
#' @param v Variance vector (length n) or NULL for constant variance.
#' @param offset Optional offset vector.
#' @param tol Convergence tolerance.
#' @param max_iter Maximum iterations.
#' @return List with components: beta, mu, eta, converged, iterations.
#' @keywords internal
.irls_loop <- function(x, y, beta_init, v = NULL, offset = NULL,
                       tol = 1e-5, max_iter = 100) {
  beta_old <- beta_init + 0.1  # perturb to enter loop
  beta_new <- beta_init
  iter <- 0

  while (max(abs(beta_new - beta_old)) > tol && iter < max_iter) {
    beta_old <- beta_new
    beta_new <- .irls_step(x, y, beta_old, v, offset)
    iter <- iter + 1
  }

  eta <- .compute_eta(x, beta_new, offset)
  mu <- .link_inv(eta)

  list(
    beta = beta_new,
    mu = as.vector(mu),
    eta = as.vector(eta),
    converged = (max(abs(beta_new - beta_old)) <= tol),
    iterations = iter
  )
}

#' Estimate Variance Function via P-Spline
#'
#' Uses mgcv::gam() with P-spline basis to estimate V(mu) nonparametrically,
#' as described in Step 2 of the estimation procedure (Equation 5).
#'
#' @param mu Fitted mean vector.
#' @param resid_sq Squared residuals (y - mu)^2.
#' @param use_log_link Whether to use log link in the GAM for variance
#'   estimation. TRUE for the standard fql, FALSE when fitting with
#'   identity link on the variance (fql_log style).
#' @return List with v (estimated variance vector) and gam_fit (the GAM object).
#' @keywords internal
.estimate_var_spline <- function(mu, resid_sq, use_log_link = TRUE) {
  dat <- data.frame(X1 = resid_sq, X2 = as.vector(mu))

  if (use_log_link) {
    # Log the squared residuals, fit with identity link, then exponentiate
    dat$X1 <- log(resid_sq)
    gam_fit <- mgcv::gam(X1 ~ s(X2, bs = "ps"),
                         family = stats::gaussian(link = "identity"),
                         data = dat)
    newd <- data.frame(X2 = as.vector(mu))
    v <- exp(stats::predict(gam_fit, newd))
  } else {
    # Fit squared residuals directly with log link
    gam_fit <- mgcv::gam(X1 ~ s(X2, bs = "ps"),
                         family = stats::gaussian(link = "log"),
                         data = dat)
    newd <- data.frame(X2 = as.vector(mu))
    v <- exp(stats::predict(gam_fit, newd))
  }

  # Ensure variance is positive
  v <- pmax(as.vector(v), .Machine$double.eps)

  list(v = v, gam_fit = gam_fit)
}

#' Estimate Variance Function via Kernel Smoother
#'
#' Uses a Gaussian kernel smoother to estimate V(mu) as an alternative
#' to P-splines when the spline method fails.
#'
#' @param mu Fitted mean vector.
#' @param resid_sq Squared residuals (y - mu)^2.
#' @param bw Bandwidth for the kernel smoother.
#' @return Estimated variance vector.
#' @keywords internal
.estimate_var_kernel <- function(mu, resid_sq, bw) {
  n <- length(mu)
  v <- numeric(n)
  for (i in seq_len(n)) {
    v[i] <- .kernel_varest(mu[i], mu, resid_sq, bw)
  }
  # Ensure variance is positive
  v <- pmax(v, .Machine$double.eps)
  v
}

#' Compute Sandwich Covariance Matrix
#'
#' Computes the sandwich (robust) covariance matrix for the FQL estimator,
#' as given in Equation (7) of the paper:
#'
#' \eqn{Var(\hat\beta) = A^{-1} B A^{-1}}
#'
#' where \eqn{A = \sum D_i D_i^T / V_i} and \eqn{B = \sum D_i D_i^T (y_i - \mu_i)^2 / V_i^2}.
#'
#' @param Di Derivative matrix (n x p).
#' @param v Variance vector (length n).
#' @param y Response vector.
#' @param mu Fitted mean vector.
#' @return Covariance matrix (p x p).
#' @keywords internal
.sandwich_vcov <- function(Di, v, y, mu) {
  p <- ncol(Di)
  n <- nrow(Di)
  v <- as.vector(v)

  # A = sum Di Di^T / v_i = crossprod(Di / sqrt(v))
  w <- 1 / sqrt(v)
  Di_w <- Di * w
  A <- crossprod(Di_w)

  # B = sum Di Di^T * (y_i - mu_i)^2 / v_i^2
  resid <- as.vector(y - mu)
  w2 <- resid / v
  Di_w2 <- Di * w2
  B <- crossprod(Di_w2)

  # Sandwich: A^{-1} B A^{-1}
  A_inv <- .safe_solve(A, diag(p), p)
  A_inv %*% B %*% A_inv
}

#' Fit FQL Model with a Specific Variance Estimation Method
#'
#' Core inner fitting routine. Alternates between estimating the variance
#' function V(mu) and re-estimating beta until convergence.
#'
#' @param x Design matrix (n x p).
#' @param y Response vector.
#' @param beta_init Initial beta estimate.
#' @param offset Optional offset vector.
#' @param var_method Variance estimation method: "spline", "kernel", or "constant".
#' @param bw Bandwidth for kernel method.
#' @param tol Convergence tolerance.
#' @param max_outer Maximum outer iterations (alternating V and beta).
#' @param max_inner Maximum inner iterations (IRLS for beta).
#' @param use_log_spline Whether to log-transform residuals for spline fitting.
#' @return List with fitting results.
#' @keywords internal
.fit_fql_core <- function(x, y, beta_init, offset = NULL,
                          var_method = "spline", bw = NULL,
                          tol = 1e-5, max_outer = 100, max_inner = 100,
                          use_log_spline = TRUE) {
  n <- nrow(x)
  p <- ncol(x)

  bn <- beta_init
  bo <- bn + 0.1
  iter_outer <- 0
  v <- rep(1, n)
  gam_fit <- NULL
  prev_change <- Inf
  damping <- 1.0  # start with no dampening

  # Accumulator for averaging when oscillation is detected
  beta_accum <- matrix(0, nrow = 0, ncol = p)

  while (iter_outer < max_outer) {
    bo <- bn

    # Compute current mu and residuals
    eta <- .compute_eta(x, bo, offset)
    mu <- .link_inv(eta)
    resid_sq <- (y - mu)^2

    # Step 2: Estimate variance function
    if (var_method == "spline") {
      var_est <- suppressWarnings(
        .estimate_var_spline(mu, resid_sq, use_log_link = use_log_spline)
      )
      v <- var_est$v
      gam_fit <- var_est$gam_fit
    } else if (var_method == "kernel") {
      v <- .estimate_var_kernel(mu, resid_sq, bw)
    } else {
      # constant: v = 1 for all subjects
      v <- rep(1, n)
    }

    # Step 3: Re-estimate beta with updated variance
    irls_result <- .irls_loop(x, y, bo, v = v, offset = offset,
                              tol = tol, max_iter = max_inner)
    bn_raw <- irls_result$beta

    # Apply dampening: blend new and old estimates to stabilize oscillation
    bn <- damping * bn_raw + (1 - damping) * bo
    iter_outer <- iter_outer + 1

    # Check convergence with relaxed tolerance for alternating optimization
    max_change <- max(abs(bn - bo))
    rel_change <- max_change / (max(abs(bo)) + tol)
    outer_tol <- max(tol * 100, 1e-4)

    if (rel_change < outer_tol) break

    # Detect oscillation: if change is not decreasing, apply dampening
    if (iter_outer > 3 && max_change > 0.8 * prev_change && damping > 0.3) {
      damping <- damping * 0.7  # increase dampening
    }
    prev_change <- max_change

    # Collect recent betas for averaging as a fallback
    if (iter_outer > max_outer * 0.5) {
      beta_accum <- rbind(beta_accum, bn)
    }
  }

  # If still oscillating after max_outer, use the average of recent betas
  if (iter_outer >= max_outer && nrow(beta_accum) > 2) {
    bn <- colMeans(beta_accum)
  }

  # Final quantities
  eta <- .compute_eta(x, bn, offset)
  mu <- as.vector(.link_inv(eta))
  resid_sq <- (y - mu)^2
  Di <- .compute_Di(mu, x)

  # Recompute A matrix for rank check
  w <- 1 / sqrt(as.vector(v))
  Di_w <- Di * w
  DV <- crossprod(Di_w)
  rank <- .matrix_rank(DV)

  # Report convergence using the same relaxed criterion
  final_rel_change <- max(abs(bn - bo)) / (max(abs(bo)) + tol)

  list(
    beta = bn,
    mu = mu,
    eta = as.vector(eta),
    v = as.vector(v),
    Di = Di,
    DV = DV,
    resid_sq = resid_sq,
    rank = rank,
    converged = (final_rel_change < outer_tol),
    iterations = iter_outer,
    gam_fit = gam_fit
  )
}
