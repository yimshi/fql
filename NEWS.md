# fql 2.0.0

This is a major rewrite of the fql package with a modern R package architecture, complete documentation, and significant improvements to code quality and usability.

## Breaking Changes

- `fql()` now returns an S3 object of class `"fql"` instead of an unnamed list. Access results via standard methods (`coef()`, `summary()`, etc.) instead of `[[1]]`, `[[2]]`.
- `fql_log()` has been removed. Use `fql(..., init.method = "zero")` instead.
- The function signature has changed: `fql(formula, data, ...)` now uses a formula interface.

## New Features

- **S3 class system**: `fql` objects work with standard R generics:
  `print()`, `summary()`, `coef()`, `vcov()`, `confint()`, `predict()`,
  `residuals()`, `fitted()`, `logLik()`, `nobs()`, `plot()`
- **`fql_multi()`**: New function for multi-taxon differential abundance analysis with FDR correction, designed for microbiome studies.
- **Diagnostic plots**: `plot.fql()` produces a 4-panel diagnostic display (residuals vs fitted, Q-Q plot, scale-location, variance function).
- **`predict()` with `newdata`**: Supports prediction on new observations.
- **Simulated example dataset**: `sim_microbiome` dataset included for demonstrations.
- **Vignette**: Comprehensive tutorial vignette ("Introduction to FQL").

## Improvements

- **Unified API**: `fql()` and `fql_log()` merged into a single `fql()` function with `init.method` parameter (`"nb"` or `"zero"`).
- **Variance method control**: New `var.method` parameter (`"auto"`, `"spline"`, `"kernel"`, `"constant"`) gives users explicit control over variance estimation.
- **Vectorized computation**: Inner loops replaced with vectorized matrix operations (`crossprod()`) for better performance.
- **Convergence stability**: Dampening strategy for the alternating optimization to handle oscillation in variance/beta estimation.
- **Formula interface**: Standard R formula with support for `offset()`.
- **Complete roxygen2 documentation**: All exported functions and methods fully documented with examples.
- **Test suite**: 27 tests covering core fitting, S3 methods, and edge cases.

## Bug Fixes

- Fixed incorrect variable reference (`b` instead of `bn`) in the constant-variance fallback path.
- Fixed missing `sqrt()` when computing standard errors from the sandwich covariance diagonal.
- Fixed p-value calculation: was using `coefficients^2/se` instead of `coefficients^2/se^2`.

## References

Shi, Y., Li, H., Wang, C., Chen, J., Jiang, H., Shih, Y.-C. T., Zhang, H., Song, Y., Feng, Y., & Liu, L. (2023). A flexible quasi-likelihood model for microbiome abundance count data. *Statistics in Medicine*, 42(25), 4632--4643. doi:10.1002/sim.9880
