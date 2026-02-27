#' fql: Flexible Quasi-Likelihood Models for Count Data
#'
#' The \pkg{fql} package fits flexible quasi-likelihood (FQL) generalized
#' linear models for count data, particularly microbiome abundance counts.
#' Unlike standard parametric approaches (Poisson, negative binomial), FQL
#' does not assume a specific distribution for the response. Instead, it
#' models the variance as an unknown smooth function of the mean, estimated
#' nonparametrically via P-splines or kernel smoothing.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{fql}}}{Fit a single FQL model.}
#'   \item{\code{\link{fql_multi}}}{Fit FQL models across multiple taxa with
#'     FDR correction.}
#' }
#'
#' @section S3 methods:
#' The following methods are available for \code{"fql"} objects:
#' \code{\link{print.fql}}, \code{\link{summary.fql}},
#' \code{\link{coef.fql}}, \code{\link{vcov.fql}},
#' \code{\link{confint.fql}}, \code{\link{predict.fql}},
#' \code{\link{residuals.fql}}, \code{\link{fitted.fql}},
#' \code{\link{plot.fql}}, \code{\link{nobs.fql}},
#' \code{\link{logLik.fql}}.
#'
#' @references
#' Shi, Y., Li, H., Wang, C., Chen, J., Jiang, H., Shih, Y.-C. T.,
#' Zhang, H., Song, Y., Feng, Y., & Liu, L. (2023).
#' A flexible quasi-likelihood model for microbiome abundance count data.
#' \emph{Statistics in Medicine}, 42(25), 4632--4643.
#' \doi{10.1002/sim.9880}
#'
#' @docType package
#' @name fql-package
#' @aliases fql-package
"_PACKAGE"
