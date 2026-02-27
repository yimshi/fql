#' Diagnostic Plots for an FQL Object
#'
#' Produces a panel of four diagnostic plots for a fitted FQL model:
#' \enumerate{
#'   \item Pearson residuals vs. fitted values
#'   \item Normal Q-Q plot of Pearson residuals
#'   \item Scale-location plot (sqrt of |Pearson residuals| vs. fitted values)
#'   \item Estimated variance function V(mu) vs. fitted values
#' }
#'
#' @param x An object of class \code{"fql"}.
#' @param which Integer vector specifying which plots to produce (1 to 4).
#'   Default is \code{1:4} (all four).
#' @param ask Logical. If \code{TRUE} and multiple plots are requested in an
#'   interactive session, the user is prompted before each new plot.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' @return Invisibly returns \code{x}.
#'
#' @details
#' \describe{
#'   \item{Plot 1 (Residuals vs Fitted)}{Shows Pearson residuals against fitted
#'     values. Should show no systematic pattern if the model is adequate.}
#'   \item{Plot 2 (Q-Q Plot)}{Normal quantile-quantile plot of Pearson
#'     residuals. Points should fall along the diagonal if the residuals are
#'     approximately normal.}
#'   \item{Plot 3 (Scale-Location)}{Shows the square root of absolute Pearson
#'     residuals against fitted values. Useful for detecting
#'     heteroscedasticity.}
#'   \item{Plot 4 (Variance Function)}{Shows the estimated variance V(mu)
#'     against fitted values, illustrating the estimated mean-variance
#'     relationship.}
#' }
#'
#' @examples
#' data(quine, package = "MASS")
#' fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
#' plot(fit)
#'
#' # Only the Q-Q plot
#' plot(fit, which = 2)
#'
#' @importFrom graphics par plot abline smoothScatter lines
#' @importFrom stats qqnorm qqline lowess residuals
#' @export
plot.fql <- function(x, which = 1:4, ask = (length(which) > 1 &&
                     grDevices::dev.interactive()), ...) {
  mu <- x$fitted.values
  pres <- residuals(x, type = "pearson")
  v <- x$variance

  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }

  if (1 %in% which) {
    graphics::plot(mu, pres,
         xlab = "Fitted values",
         ylab = "Pearson residuals",
         main = "Residuals vs Fitted",
         pch = 1, col = grDevices::adjustcolor("black", alpha.f = 0.5), ...)
    graphics::abline(h = 0, lty = 2, col = "grey50")
    lo <- stats::lowess(mu, pres)
    graphics::lines(lo, col = "red", lwd = 1.5)
  }

  if (2 %in% which) {
    stats::qqnorm(pres, main = "Q-Q Plot of Pearson Residuals",
                  pch = 1, col = grDevices::adjustcolor("black", alpha.f = 0.5), ...)
    stats::qqline(pres, col = "red", lwd = 1.5)
  }

  if (3 %in% which) {
    sqrt_abs_pres <- sqrt(abs(pres))
    graphics::plot(mu, sqrt_abs_pres,
         xlab = "Fitted values",
         ylab = expression(sqrt("|Pearson residuals|")),
         main = "Scale-Location",
         pch = 1, col = grDevices::adjustcolor("black", alpha.f = 0.5), ...)
    lo <- stats::lowess(mu, sqrt_abs_pres)
    graphics::lines(lo, col = "red", lwd = 1.5)
  }

  if (4 %in% which) {
    ord <- order(mu)
    graphics::plot(mu[ord], v[ord],
         xlab = expression(hat(mu)),
         ylab = expression(hat(V)(mu)),
         main = "Estimated Variance Function",
         type = "p", pch = 1,
         col = grDevices::adjustcolor("steelblue", alpha.f = 0.6), ...)
    lo <- stats::lowess(mu[ord], v[ord])
    graphics::lines(lo, col = "red", lwd = 2)
  }

  invisible(x)
}
