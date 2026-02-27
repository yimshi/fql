# fql

`fql` provides Flexible Quasi-Likelihood (FQL) regression tools for microbiome and
other count data settings with heteroskedasticity.

## Installation

```r
# install.packages("remotes")
remotes::install_github("yimshi/fql")
```

## Quick start

```r
library(fql)
library(MASS)

data("quine", package = "MASS")
fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine, na.action = na.omit)

fit$estimation
fit$covariance
fit$meanerror
fit$step
```

## Release resources

- Changelog: `NEWS.md`
- Citation information: `citation("fql")`
- Issue tracker: <https://github.com/yimshi/fql/issues>
