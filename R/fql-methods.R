#' Print an FQL Object
#'
#' Prints a brief summary of a fitted FQL model, including the call and
#' estimated coefficients.
#'
#' @param x An object of class \code{"fql"}.
#' @param digits Number of significant digits. Default is
#'   \code{max(3, getOption("digits") - 3)}.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' print(fit)
#'
#' @export
print.fql <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nFlexible Quasi-Likelihood Model (FQL)\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(round(x$coefficients, digits))
  cat("\nVariance estimation method:", x$var.method.used, "\n")
  cat("Number of observations:", x$n, "\n")
  if (!x$converged) {
    cat("WARNING: Algorithm did not converge.\n")
  }
  invisible(x)
}

#' Summarize an FQL Object
#'
#' Produces a detailed summary of a fitted FQL model, including coefficient
#' estimates, standard errors, Wald z-statistics, and p-values with
#' significance stars.
#'
#' @param object An object of class \code{"fql"}.
#' @param ... Additional arguments (currently unused).
#' @return An object of class \code{"summary.fql"}, which is a list containing:
#'   \describe{
#'     \item{call}{The model call.}
#'     \item{coefficients}{A matrix with columns for Estimate, Std. Error,
#'       z value, and Pr(>|z|).}
#'     \item{var.method.used}{Variance estimation method used.}
#'     \item{n}{Number of observations.}
#'     \item{df.residual}{Residual degrees of freedom.}
#'     \item{converged}{Whether the algorithm converged.}
#'     \item{iterations}{Number of iterations.}
#'     \item{link}{Link function used.}
#'     \item{mean.residual}{Mean of squared residuals.}
#'   }
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' summary(fit)
#'
#' @export
summary.fql <- function(object, ...) {
  z <- object$coefficients / object$se
  coef_table <- cbind(
    Estimate = object$coefficients,
    `Std. Error` = object$se,
    `z value` = z,
    `Pr(>|z|)` = object$p.value
  )
  rownames(coef_table) <- names(object$coefficients)

  out <- list(
    call = object$call,
    coefficients = coef_table,
    var.method.used = object$var.method.used,
    n = object$n,
    df.residual = object$df.residual,
    converged = object$converged,
    iterations = object$iterations,
    link = object$link,
    mean.residual = mean(object$residuals^2)
  )
  class(out) <- "summary.fql"
  out
}

#' Print Summary of an FQL Object
#'
#' @param x An object of class \code{"summary.fql"}.
#' @param digits Number of significant digits.
#' @param signif.stars Logical: print significance stars?
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns \code{x}.
#' @export
print.summary.fql <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars", TRUE),
                               ...) {
  cat("\nFlexible Quasi-Likelihood Model (FQL)\n\n")
  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars,
                      na.print = "NA", ...)

  cat("\n---\n")
  cat("Link function:", x$link, "\n")
  cat("Variance estimation:", x$var.method.used, "\n")
  cat("Number of observations:", x$n, "\n")
  cat("Residual degrees of freedom:", x$df.residual, "\n")
  cat("Mean squared residual:", format(x$mean.residual, digits = digits), "\n")
  if (!x$converged) {
    cat("\nWARNING: Algorithm did not converge after ", x$iterations,
        " iterations.\n", sep = "")
  }
  invisible(x)
}

#' Extract Coefficients from an FQL Object
#'
#' @param object An object of class \code{"fql"}.
#' @param ... Additional arguments (currently unused).
#' @return Named vector of estimated coefficients.
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' coef(fit)
#'
#' @export
coef.fql <- function(object, ...) {
  object$coefficients
}

#' Extract Variance-Covariance Matrix from an FQL Object
#'
#' Returns the sandwich variance-covariance matrix of the estimated
#' regression coefficients.
#'
#' @param object An object of class \code{"fql"}.
#' @param ... Additional arguments (currently unused).
#' @return The variance-covariance matrix.
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' vcov(fit)
#'
#' @export
vcov.fql <- function(object, ...) {
  object$vcov
}

#' Confidence Intervals for FQL Coefficients
#'
#' Computes Wald-type confidence intervals for the regression coefficients.
#'
#' @param object An object of class \code{"fql"}.
#' @param parm A specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param level The confidence level. Default is 0.95.
#' @param ... Additional arguments (currently unused).
#' @return A matrix with columns for the lower and upper confidence limits.
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' confint(fit)
#' confint(fit, level = 0.99)
#'
#' @export
confint.fql <- function(object, parm, level = 0.95, ...) {
  cf <- object$coefficients
  se <- object$se
  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)

  if (missing(parm)) {
    parm <- names(cf)
  } else if (is.numeric(parm)) {
    parm <- names(cf)[parm]
  }

  ci <- cbind(cf[parm] - z * se[parm], cf[parm] + z * se[parm])
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
  ci
}

#' Predict from an FQL Object
#'
#' Obtains predictions from a fitted FQL model, optionally for new data.
#'
#' @param object An object of class \code{"fql"}.
#' @param newdata Optional data frame for prediction. If omitted, predictions
#'   are for the original data.
#' @param type The type of prediction: \code{"response"} (default) returns
#'   predictions on the response scale (mu), \code{"link"} returns predictions
#'   on the linear predictor scale (eta).
#' @param ... Additional arguments (currently unused).
#' @return A vector of predictions.
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' head(predict(fit, type = "response"))
#' head(predict(fit, type = "link"))
#'
#' @export
predict.fql <- function(object, newdata = NULL, type = c("response", "link"),
                        ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    # Build design matrix for new data
    # Use the stored design matrix column names to ensure consistency
    tt <- stats::delete.response(stats::terms(object$formula))
    mf <- stats::model.frame(tt, data = newdata)
    x_new <- stats::model.matrix(tt, mf)

    # Handle offset
    offset_new <- as.vector(stats::model.offset(mf))

    eta <- .compute_eta(x_new, object$coefficients, offset_new)
  }

  if (type == "response") {
    .link_inv(eta)
  } else {
    eta
  }
}

#' Extract Residuals from an FQL Object
#'
#' @param object An object of class \code{"fql"}.
#' @param type The type of residuals: \code{"response"} (default, y - mu),
#'   \code{"pearson"} ((y - mu) / sqrt(V(mu))), or \code{"working"}
#'   ((y - mu) / mu for log link).
#' @param ... Additional arguments (currently unused).
#' @return A vector of residuals.
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' head(residuals(fit, type = "pearson"))
#'
#' @export
residuals.fql <- function(object, type = c("response", "pearson", "working"),
                          ...) {
  type <- match.arg(type)

  switch(type,
    response = object$residuals,
    pearson = object$residuals / sqrt(object$variance),
    working = object$residuals / object$fitted.values
  )
}

#' Extract Fitted Values from an FQL Object
#'
#' @param object An object of class \code{"fql"}.
#' @param ... Additional arguments (currently unused).
#' @return A vector of fitted mean values.
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' head(fitted(fit))
#'
#' @export
fitted.fql <- function(object, ...) {
  object$fitted.values
}

#' Extract the Number of Observations from an FQL Object
#'
#' @param object An object of class \code{"fql"}.
#' @param ... Additional arguments (currently unused).
#' @return Integer: the number of observations.
#' @importFrom stats nobs
#' @export
nobs.fql <- function(object, ...) {
  object$n
}

#' Compute Quasi-Log-Likelihood for an FQL Object
#'
#' Computes the quasi-log-likelihood value, defined as
#' \eqn{Q^*(\mu, y) = \sum_i \int_{y_i}^{\mu_i} (y_i - t) / V(t) dt}.
#' This is approximated numerically.
#'
#' @param object An object of class \code{"fql"}.
#' @param ... Additional arguments (currently unused).
#' @return An object of class \code{"logLik"} with the quasi-log-likelihood value.
#' @importFrom stats logLik
#' @export
logLik.fql <- function(object, ...) {
  # Approximate the quasi-log-likelihood
  # Q* = sum_i integral_{y_i}^{mu_i} (y_i - t) / V(t) dt
  # For simplicity, approximate as: -0.5 * sum((y - mu)^2 / V)
  y <- object$y
  mu <- object$fitted.values
  v <- object$variance

  val <- -0.5 * sum((y - mu)^2 / v)
  attr(val, "df") <- length(object$coefficients)
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"
  val
}
