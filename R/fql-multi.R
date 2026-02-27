#' Fit FQL Models Across Multiple Taxa
#'
#' Applies the FQL model to each column (taxon) of an OTU/abundance count
#' table, then adjusts p-values for multiple testing. This implements the
#' analysis workflow described in Section 4 of Shi et al. (2023).
#'
#' @param otu_table A matrix or data frame of abundance counts, where rows
#'   are samples and columns are taxa (OTUs). Column names will be used as
#'   taxon identifiers in the output.
#' @param formula A one-sided formula specifying the covariates (without the
#'   response), e.g., \code{~ group + age + offset(log_total)}. The response
#'   is automatically set to each column of \code{otu_table}.
#' @param data A data frame containing the covariates referenced in
#'   \code{formula}. Must have the same number of rows as \code{otu_table}.
#' @param p.adjust.method Method for p-value adjustment. Any method accepted
#'   by \code{\link[stats]{p.adjust}} can be used. Default is \code{"BH"}
#'   (Benjamini-Hochberg for FDR control).
#' @param prevalence.filter Minimum prevalence threshold (proportion of
#'   non-zero samples). Taxa with prevalence below this value are excluded.
#'   Default is 0.25 (25\%).
#' @param verbose Logical. If \code{TRUE}, prints progress. Default is \code{TRUE}.
#' @param ... Additional arguments passed to \code{\link{fql}}.
#'
#' @return A data frame with one row per taxon and columns:
#'   \describe{
#'     \item{taxon}{Taxon name (from column names of otu_table).}
#'     \item{prevalence}{Proportion of non-zero samples.}
#'     \item{mean_abundance}{Mean abundance across all samples.}
#'     \item{covariate}{Name of each covariate (one row per taxon-covariate
#'       combination, excluding intercept).}
#'     \item{estimate}{Estimated coefficient.}
#'     \item{se}{Standard error.}
#'     \item{z_value}{Wald z-statistic.}
#'     \item{p_value}{Raw p-value.}
#'     \item{p_adjusted}{Adjusted p-value.}
#'     \item{var_method}{Variance estimation method used.}
#'     \item{converged}{Whether the model converged.}
#'   }
#'   The data frame is sorted by adjusted p-value.
#'
#' @details
#' For each taxon, the function constructs a two-sided formula with the
#' taxon count as the response and the covariates from \code{formula},
#' then calls \code{\link{fql}} to fit the model. Taxa failing to converge
#' or producing errors are reported with \code{NA} values.
#'
#' The prevalence filter removes taxa with too many zeros before fitting,
#' as extremely sparse taxa may not yield meaningful results with
#' one-part models (see Discussion in Shi et al., 2023).
#'
#' @examples
#' \donttest{
#' # Simulate a small OTU table
#' set.seed(42)
#' n <- 100
#' otu <- matrix(rnbinom(n * 5, mu = 10, size = 2), nrow = n)
#' colnames(otu) <- paste0("OTU", 1:5)
#' covariates <- data.frame(group = factor(rep(c("A", "B"), each = n/2)))
#'
#' results <- fql_multi(otu, ~ group, data = covariates,
#'                      prevalence.filter = 0)
#' print(results)
#' }
#'
#' @references
#' Shi, Y., Li, H., Wang, C., Chen, J., Jiang, H., Shih, Y.-C. T.,
#' Zhang, H., Song, Y., Feng, Y., & Liu, L. (2023).
#' A flexible quasi-likelihood model for microbiome abundance count data.
#' \emph{Statistics in Medicine}, 42(25), 4632--4643.
#' \doi{10.1002/sim.9880}
#'
#' @seealso \code{\link{fql}}, \code{\link[stats]{p.adjust}}
#'
#' @importFrom stats p.adjust
#' @export
fql_multi <- function(otu_table, formula, data,
                      p.adjust.method = "BH",
                      prevalence.filter = 0.25,
                      verbose = TRUE, ...) {

  # Input validation
  if (!is.matrix(otu_table) && !is.data.frame(otu_table)) {
    stop("'otu_table' must be a matrix or data frame.", call. = FALSE)
  }
  otu_table <- as.matrix(otu_table)

  if (nrow(otu_table) != nrow(data)) {
    stop("'otu_table' and 'data' must have the same number of rows.",
         call. = FALSE)
  }

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a one-sided formula (e.g., ~ group + age).",
         call. = FALSE)
  }

  # Ensure formula is one-sided
  if (length(formula) == 3) {
    warning("Two-sided formula detected. Ignoring the left-hand side.",
            call. = FALSE)
    formula <- formula[-2]
  }

  # Taxon names
  taxa_names <- colnames(otu_table)
  if (is.null(taxa_names)) {
    taxa_names <- paste0("Taxon", seq_len(ncol(otu_table)))
  }

  # Prevalence filter
  prevalence <- colMeans(otu_table > 0)
  keep <- prevalence >= prevalence.filter
  n_filtered <- sum(!keep)

  if (n_filtered > 0 && verbose) {
    message(sprintf("Filtered %d taxa with prevalence < %.0f%% (%d taxa remaining).",
                    n_filtered, prevalence.filter * 100, sum(keep)))
  }

  if (sum(keep) == 0) {
    stop("No taxa passed the prevalence filter.", call. = FALSE)
  }

  otu_filtered <- otu_table[, keep, drop = FALSE]
  taxa_filtered <- taxa_names[keep]
  prevalence_filtered <- prevalence[keep]

  # Fit FQL to each taxon
  n_taxa <- ncol(otu_filtered)
  results_list <- vector("list", n_taxa)

  for (j in seq_len(n_taxa)) {
    taxon_name <- taxa_filtered[j]
    if (verbose) {
      message(sprintf("[%d/%d] Fitting %s...", j, n_taxa, taxon_name))
    }

    # Create response column in data
    fit_data <- data
    fit_data$.response <- otu_filtered[, j]

    # Build two-sided formula
    rhs <- as.character(formula)[2]
    fit_formula <- stats::as.formula(paste(".response ~", rhs))

    # Fit model
    result <- tryCatch({
      fit <- fql(fit_formula, data = fit_data, verbose = FALSE, ...)

      # Extract results for non-intercept coefficients
      coef_names <- names(fit$coefficients)
      non_intercept <- coef_names[coef_names != "(Intercept)"]

      data.frame(
        taxon = taxon_name,
        prevalence = prevalence_filtered[j],
        mean_abundance = mean(otu_filtered[, j]),
        covariate = non_intercept,
        estimate = fit$coefficients[non_intercept],
        se = fit$se[non_intercept],
        z_value = fit$coefficients[non_intercept] / fit$se[non_intercept],
        p_value = fit$p.value[non_intercept],
        var_method = fit$var.method.used,
        converged = fit$converged,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      if (verbose) {
        message(sprintf("  Warning: failed for %s (%s)", taxon_name,
                        conditionMessage(e)))
      }
      # Return NA row for failed taxa
      coef_names <- tryCatch({
        mf <- stats::model.frame(formula, data = data)
        mt <- attr(mf, "terms")
        x <- stats::model.matrix(mt, mf)
        cn <- colnames(x)
        cn[cn != "(Intercept)"]
      }, error = function(e2) "unknown")

      data.frame(
        taxon = taxon_name,
        prevalence = prevalence_filtered[j],
        mean_abundance = mean(otu_filtered[, j]),
        covariate = coef_names,
        estimate = NA_real_,
        se = NA_real_,
        z_value = NA_real_,
        p_value = NA_real_,
        var_method = NA_character_,
        converged = FALSE,
        row.names = NULL,
        stringsAsFactors = FALSE
      )
    })

    results_list[[j]] <- result
  }

  # Combine results
  results <- do.call(rbind, results_list)

  # Adjust p-values for multiple testing (within each covariate)
  covariates_unique <- unique(results$covariate)
  results$p_adjusted <- NA_real_

  for (cov in covariates_unique) {
    idx <- results$covariate == cov
    results$p_adjusted[idx] <- stats::p.adjust(results$p_value[idx],
                                               method = p.adjust.method)
  }

  # Sort by adjusted p-value
  results <- results[order(results$p_adjusted, na.last = TRUE), ]
  rownames(results) <- NULL

  results
}
