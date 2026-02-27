#' Construct an FQL Object
#'
#' Internal constructor for creating objects of class \code{"fql"}.
#'
#' @param coefficients Named vector of regression coefficients.
#' @param se Named vector of standard errors.
#' @param p.value Named vector of p-values.
#' @param vcov Variance-covariance matrix.
#' @param fitted.values Fitted mean values.
#' @param residuals Response residuals.
#' @param variance Estimated variance function values.
#' @param var.method.used Character describing variance method used.
#' @param linear.predictors Linear predictor values.
#' @param converged Logical: did the algorithm converge?
#' @param iterations Number of outer iterations.
#' @param df.residual Residual degrees of freedom.
#' @param n Number of observations.
#' @param formula Model formula.
#' @param call Matched call.
#' @param link Link function name.
#' @param x Design matrix.
#' @param y Response vector.
#' @param offset Offset vector or NULL.
#' @return An object of class \code{"fql"}.
#' @keywords internal
.new_fql <- function(coefficients, se, p.value, vcov, fitted.values,
                     residuals, variance, var.method.used, linear.predictors,
                     converged, iterations, df.residual, n, formula, call,
                     link, x, y, offset) {
  structure(
    list(
      coefficients = coefficients,
      se = se,
      p.value = p.value,
      vcov = vcov,
      fitted.values = fitted.values,
      residuals = residuals,
      variance = variance,
      var.method.used = var.method.used,
      linear.predictors = linear.predictors,
      converged = converged,
      iterations = iterations,
      df.residual = df.residual,
      n = n,
      formula = formula,
      call = call,
      link = link,
      x = x,
      y = y,
      offset = offset
    ),
    class = "fql"
  )
}
