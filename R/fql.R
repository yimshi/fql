#' Fit a Flexible Quasi-Likelihood Model
#'
#' Fits a flexible quasi-likelihood (FQL) generalized linear model for count
#' data with heteroscedastic variance. The variance is modeled as an unknown
#' smooth function of the mean, estimated nonparametrically via P-splines or
#' kernel smoothing. This avoids distributional assumptions while properly
#' accounting for overdispersion.
#'
#' @param formula An object of class \code{\link[stats]{formula}}: a symbolic
#'   description of the model to be fitted. The response should be a
#'   non-negative count variable.
#' @param data A data frame containing the variables in the model.
#' @param init.method Method for initializing beta coefficients. \code{"nb"}
#'   (default) initializes from a negative binomial GLM fit via
#'   \code{\link[MASS]{glm.nb}}. \code{"zero"} initializes all coefficients
#'   to zero.
#' @param var.method Method for estimating the variance function V(mu).
#'   \code{"auto"} (default) tries methods in order: spline, kernel (bw+1),
#'   kernel (bw+3), constant, using the first that succeeds. \code{"spline"}
#'   uses P-spline via \code{\link[mgcv]{gam}}. \code{"kernel"} uses Gaussian
#'   kernel smoothing. \code{"constant"} sets V(mu) = 1 (equivalent to
#'   Poisson-like working variance).
#' @param tol Convergence tolerance for the iterative algorithm. Default is
#'   \code{1e-5}.
#' @param max.iter Maximum number of outer iterations (alternating between
#'   variance estimation and beta estimation). Default is \code{100}.
#' @param verbose Logical. If \code{TRUE}, prints progress messages during
#'   fitting. Default is \code{FALSE}.
#'
#' @return An object of class \code{"fql"}, which is a list containing:
#'   \describe{
#'     \item{coefficients}{Named vector of estimated regression coefficients.}
#'     \item{se}{Named vector of sandwich standard errors.}
#'     \item{p.value}{Named vector of Wald test p-values.}
#'     \item{vcov}{Sandwich variance-covariance matrix.}
#'     \item{fitted.values}{Estimated means (mu).}
#'     \item{residuals}{Response residuals (y - mu).}
#'     \item{variance}{Estimated variance function values V(mu_i).}
#'     \item{var.method.used}{Character string indicating which variance
#'       estimation method was used: \code{"spline"}, \code{"kernel"},
#'       or \code{"constant"}.}
#'     \item{linear.predictors}{Linear predictor values (eta = X beta).}
#'     \item{converged}{Logical indicating whether the algorithm converged.}
#'     \item{iterations}{Number of outer iterations performed.}
#'     \item{df.residual}{Residual degrees of freedom.}
#'     \item{n}{Number of observations.}
#'     \item{formula}{The model formula.}
#'     \item{call}{The matched call.}
#'     \item{link}{The link function used (currently always \code{"log"}).}
#'     \item{x}{The design matrix.}
#'     \item{y}{The response vector.}
#'   }
#'
#' @details
#' The FQL model assumes:
#' \deqn{\log(\mu_i) = X_i^T \beta}
#' \deqn{Var(Y_i) = V(\mu_i)}
#' where \eqn{V(\cdot)} is an unknown smooth function estimated nonparametrically.
#'
#' The estimation alternates between:
#' \enumerate{
#'   \item Estimating \eqn{V(\mu)} by fitting a P-spline (or kernel smoother)
#'     to the squared residuals as a function of \eqn{\mu}.
#'   \item Re-estimating \eqn{\beta} via IRLS with the estimated variance weights.
#' }
#'
#' Inference uses a sandwich covariance estimator that is robust to
#' misspecification of the variance function.
#'
#' When \code{var.method = "auto"}, the function tries variance estimation
#' methods in the following order, using the first that succeeds:
#' \enumerate{
#'   \item P-spline via \code{mgcv::gam()}
#'   \item Kernel smoother with bandwidth \code{bw + 1}
#'   \item Kernel smoother with bandwidth \code{bw + 3}
#'   \item Constant variance (V = 1)
#' }
#' This fallback strategy handles numerical difficulties that can arise with
#' sparse or extreme data.
#'
#' @references
#' Shi, Y., Li, H., Wang, C., Chen, J., Jiang, H., Shih, Y.-C. T.,
#' Zhang, H., Song, Y., Feng, Y., & Liu, L. (2023).
#' A flexible quasi-likelihood model for microbiome abundance count data.
#' \emph{Statistics in Medicine}, 42(25), 4632--4643.
#' \doi{10.1002/sim.9880}
#'
#' @examples
#' # Fit FQL model to the quine dataset
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' summary(fit)
#'
#' # Compare with negative binomial GLM
#' nb_fit <- MASS::glm.nb(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' cbind(FQL = coef(fit), NB = coef(nb_fit))
#'
#' # Predictions
#' head(predict(fit, type = "response"))
#'
#' # Confidence intervals
#' confint(fit)
#'
#' @seealso \code{\link{summary.fql}}, \code{\link{predict.fql}},
#'   \code{\link{plot.fql}}, \code{\link{fql_multi}}
#'
#' @importFrom MASS glm.nb
#' @importFrom mgcv gam s
#' @importFrom stats model.frame model.response model.matrix model.offset
#'   gaussian predict pchisq qnorm
#' @export
fql <- function(formula, data,
                init.method = c("nb", "zero"),
                var.method = c("auto", "spline", "kernel", "constant"),
                tol = 1e-5, max.iter = 100, verbose = FALSE) {

  # Capture the call
  cl <- match.call()

  # Match arguments
  init.method <- match.arg(init.method)
  var.method <- match.arg(var.method)

 # --- Input validation ---
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (tol <= 0) {
    stop("'tol' must be a positive number.", call. = FALSE)
  }
  if (max.iter < 1) {
    stop("'max.iter' must be at least 1.", call. = FALSE)
  }

  # --- Parse formula ---
  mf <- stats::model.frame(formula = formula, data = data)
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf, "numeric")
  x <- stats::model.matrix(mt, mf)
  offset <- as.vector(stats::model.offset(mf))

  n <- nrow(x)
  p <- ncol(x)

  if (!is.null(offset) && length(offset) != n) {
    stop(sprintf("Number of offsets (%d) must equal number of observations (%d).",
                 length(offset), n), call. = FALSE)
  }

  if (any(y < 0)) {
    warning("Response variable contains negative values. ",
            "FQL is designed for non-negative count data.", call. = FALSE)
  }

  # Bandwidth for kernel methods (Silverman's rule of thumb)
  bw <- (4 / 3)^(1 / 5) * n^(-1 / 5)

  # --- Step 1: Initialize beta ---
  if (verbose) message("Step 1: Initializing beta (method: ", init.method, ")")

  if (init.method == "nb") {
    nb_fit <- tryCatch(
      MASS::glm.nb(formula = formula, link = "log", data = data),
      error = function(e) NULL
    )
    if (!is.null(nb_fit)) {
      beta_init <- stats::coef(nb_fit)
    } else {
      if (verbose) message("  NB initialization failed, falling back to zeros.")
      beta_init <- rep(0, p)
    }
  } else {
    beta_init <- rep(0, p)
  }

  # Initial IRLS with constant variance to get starting beta
  if (verbose) message("  Running initial IRLS with constant variance...")
  init_result <- .irls_loop(x, y, beta_init, v = NULL, offset = offset,
                            tol = tol, max_iter = max.iter)
  beta_start <- init_result$beta

  # --- Steps 2-3: Iterate variance estimation and beta estimation ---
  fit <- NULL
  var_method_used <- NULL

  if (var.method == "auto") {
    # Try methods in order: spline -> kernel(bw+1) -> kernel(bw+3) -> constant
    methods_to_try <- list(
      list(method = "spline", bw = NULL, label = "spline"),
      list(method = "kernel", bw = bw + 1, label = "kernel (bw+1)"),
      list(method = "kernel", bw = bw + 3, label = "kernel (bw+3)"),
      list(method = "constant", bw = NULL, label = "constant")
    )

    for (m in methods_to_try) {
      if (verbose) message("  Trying variance method: ", m$label)

      result <- tryCatch({
        f <- .fit_fql_core(x, y, beta_start, offset = offset,
                          var_method = m$method, bw = m$bw,
                          tol = tol, max_outer = max.iter, max_inner = max.iter,
                          use_log_spline = TRUE)
        # Check if the DV matrix is full rank
        if (f$rank < p) {
          if (verbose) message("    Rank deficient (", f$rank, " < ", p, "), trying next method.")
          NULL
        } else {
          f
        }
      }, error = function(e) {
        if (verbose) message("    Failed: ", conditionMessage(e))
        NULL
      })

      if (!is.null(result)) {
        fit <- result
        var_method_used <- m$label
        if (verbose) message("  Success with: ", var_method_used)
        break
      }
    }

    if (is.null(fit)) {
      stop("All variance estimation methods failed. ",
           "Consider checking your data for extreme values or collinearity.",
           call. = FALSE)
    }
  } else {
    # Use specified method
    bw_use <- switch(var.method,
      kernel = bw + 1,
      NULL
    )

    if (verbose) message("  Using variance method: ", var.method)
    fit <- .fit_fql_core(x, y, beta_start, offset = offset,
                        var_method = var.method, bw = bw_use,
                        tol = tol, max_outer = max.iter, max_inner = max.iter,
                        use_log_spline = TRUE)
    var_method_used <- var.method
  }

  # --- Compute inference ---
  Di <- .compute_Di(fit$mu, x)
  vcov_mat <- .sandwich_vcov(Di, fit$v, y, fit$mu)

  beta <- fit$beta
  se <- sqrt(diag(vcov_mat))
  z_stat <- beta / se
  p_value <- 1 - stats::pchisq(z_stat^2, df = 1)

  # Name the results
  coef_names <- colnames(x)
  names(beta) <- coef_names
  names(se) <- coef_names
  names(p_value) <- coef_names
  rownames(vcov_mat) <- coef_names
  colnames(vcov_mat) <- coef_names

  # --- Build and return the S3 object ---
  .new_fql(
    coefficients = beta,
    se = se,
    p.value = p_value,
    vcov = vcov_mat,
    fitted.values = fit$mu,
    residuals = as.vector(y - fit$mu),
    variance = fit$v,
    var.method.used = var_method_used,
    linear.predictors = fit$eta,
    converged = fit$converged,
    iterations = fit$iterations,
    df.residual = n - p,
    n = n,
    formula = formula,
    call = cl,
    link = "log",
    x = x,
    y = as.vector(y),
    offset = offset
  )
}
