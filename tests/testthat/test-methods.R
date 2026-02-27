test_that("coef.fql returns named vector", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  cf <- coef(fit)
  expect_true(is.numeric(cf))
  expect_true(!is.null(names(cf)))
  expect_length(cf, 3)
})

test_that("vcov.fql returns correct dimension matrix", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  v <- vcov(fit)
  expect_true(is.matrix(v))
  expect_equal(dim(v), c(3, 3))
  expect_true(!is.null(rownames(v)))
  # Should be positive definite (all positive diagonal)
  expect_true(all(diag(v) > 0))
})

test_that("confint.fql returns correct structure", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  ci <- confint(fit)
  expect_true(is.matrix(ci))
  expect_equal(nrow(ci), 3)
  expect_equal(ncol(ci), 2)
  # Lower < Upper
  expect_true(all(ci[, 1] < ci[, 2]))

  # Custom level
  ci99 <- confint(fit, level = 0.99)
  # 99% CI should be wider
  expect_true(all((ci99[, 2] - ci99[, 1]) >= (ci[, 2] - ci[, 1])))
})

test_that("predict.fql works with and without newdata", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  # Without newdata
  pred_resp <- predict(fit, type = "response")
  pred_link <- predict(fit, type = "link")

  expect_length(pred_resp, nrow(quine))
  expect_length(pred_link, nrow(quine))
  expect_true(all(pred_resp > 0))  # log link -> positive predictions
  expect_equal(pred_resp, exp(pred_link), tolerance = 1e-10)

  # With newdata
  newdata <- quine[1:5, ]
  pred_new <- predict(fit, newdata = newdata, type = "response")
  expect_length(pred_new, 5)
})

test_that("residuals.fql returns correct types", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  r_resp <- residuals(fit, type = "response")
  r_pear <- residuals(fit, type = "pearson")
  r_work <- residuals(fit, type = "working")

  expect_length(r_resp, nrow(quine))
  expect_length(r_pear, nrow(quine))
  expect_length(r_work, nrow(quine))

  # Response residuals = y - mu
  expect_equal(r_resp, fit$y - fit$fitted.values)
})

test_that("fitted.fql returns fitted values", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  fv <- fitted(fit)
  expect_length(fv, nrow(quine))
  expect_true(all(fv > 0))
  expect_equal(fv, fit$fitted.values)
})

test_that("nobs.fql returns correct count", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  expect_equal(nobs(fit), nrow(quine))
})

test_that("logLik.fql returns logLik object", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_true(is.finite(as.numeric(ll)))
})

test_that("print.fql runs without error", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  expect_output(print(fit), "Flexible Quasi-Likelihood")
  expect_output(print(fit), "Coefficients")
})

test_that("summary.fql produces correct output", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  s <- summary(fit)
  expect_s3_class(s, "summary.fql")
  expect_true(is.matrix(s$coefficients))
  expect_equal(ncol(s$coefficients), 4)
  expect_equal(colnames(s$coefficients),
               c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))

  # print.summary.fql
  expect_output(print(s), "Flexible Quasi-Likelihood")
  expect_output(print(s), "Signif")
})

test_that("plot.fql runs without error", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine)

  # Plotting should not error
  pdf(NULL)  # null device to suppress actual plot output
  expect_invisible(plot(fit, which = 1:4, ask = FALSE))
  dev.off()
})
