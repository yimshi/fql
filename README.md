# fql: Flexible Quasi-Likelihood Models for Count Data

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**fql** fits flexible quasi-likelihood (FQL) generalized linear models for count
data with heteroscedastic variance. The variance is modeled as an unknown smooth
function of the mean, estimated nonparametrically via P-splines or kernel
smoothing. This is particularly useful for microbiome abundance count data.

## Installation

```r
# Install from GitHub
devtools::install_github("yimshi/fql")
```

## Quick Start

```r
library(fql)
data(quine, package = "MASS")

# Fit FQL model
fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)
summary(fit)

# S3 methods work as expected
coef(fit)
confint(fit)
predict(fit, type = "response")
plot(fit)
```

## Multi-Taxon Analysis

```r
data(sim_microbiome)

results <- fql_multi(
  otu_table = sim_microbiome$otu_table,
  formula = ~ group + age + offset(log(total_reads)),
  data = sim_microbiome$sample_data,
  p.adjust.method = "BH"
)
```

## Citation

Shi, Y., Li, H., Wang, C., Chen, J., Jiang, H., Shih, Y.-C. T., Zhang, H.,
Song, Y., Feng, Y., & Liu, L. (2023). A flexible quasi-likelihood model for
microbiome abundance count data. *Statistics in Medicine*, 42(25), 4632-4643.
[doi:10.1002/sim.9880](https://doi.org/10.1002/sim.9880)
