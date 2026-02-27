#' Simulated Microbiome Abundance Data
#'
#' A simulated dataset demonstrating microbiome abundance count data with
#' known differential abundance patterns. Generated from negative binomial
#' distributions with varying dispersion parameters and effect sizes.
#'
#' @format A list with three components:
#' \describe{
#'   \item{otu_table}{A 200 x 10 integer matrix of abundance counts. Rows are
#'     samples, columns are OTUs (OTU1 through OTU10). OTUs 1-3 have true
#'     differential abundance between groups.}
#'   \item{sample_data}{A data frame with 200 rows and 3 columns:
#'     \describe{
#'       \item{group}{Factor with levels "Control" and "Treatment" (100 each).}
#'       \item{age}{Numeric, simulated ages (mean 50, sd 10).}
#'       \item{total_reads}{Integer, total read count per sample (sum across
#'         all OTUs).}
#'     }
#'   }
#'   \item{taxa_info}{A data frame with 10 rows describing each OTU:
#'     \describe{
#'       \item{taxon}{Character, OTU name.}
#'       \item{differential}{Logical, whether the OTU has true differential
#'         abundance.}
#'       \item{true_effect}{Numeric, the true group effect (log scale).
#'         Non-zero for OTUs 1-3.}
#'     }
#'   }
#' }
#'
#' @source Simulated using negative binomial distributions. See
#'   \code{data-raw/simulate-data.R} for the generation script.
#'
#' @examples
#' data(sim_microbiome)
#'
#' # Examine the structure
#' str(sim_microbiome)
#'
#' # Fit FQL to a single OTU
#' dat <- sim_microbiome$sample_data
#' dat$count <- sim_microbiome$otu_table[, "OTU1"]
#' fit <- fql(count ~ group + age, data = dat)
#' summary(fit)
#'
"sim_microbiome"
