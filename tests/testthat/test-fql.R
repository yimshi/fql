test_that("fql fits basic model with quine data", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex + Age + Lrn, data = quine)

  expect_s3_class(fit, "fql")
  expect_length(fit$coefficients, 7)  # intercept + 6 dummy vars
  expect_equal(fit$n, nrow(quine))
  expect_true(fit$converged)
  expect_equal(fit$link, "log")
})

test_that("fql with init.method = 'zero' works", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine, init.method = "zero")

  expect_s3_class(fit, "fql")
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("fql with var.method = 'constant' works", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex, data = quine, var.method = "constant")

  expect_s3_class(fit, "fql")
  expect_true(all(fit$variance == 1))
})

test_that("fql handles single covariate", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth, data = quine)

  expect_s3_class(fit, "fql")
  expect_length(fit$coefficients, 2)
})

test_that("fql handles offset", {
  data(quine, package = "MASS")
  quine$log_offset <- log(quine$Days + 1)
  fit <- fql(Days ~ Eth + Sex + offset(log_offset), data = quine)

  expect_s3_class(fit, "fql")
  expect_false(is.null(fit$offset))
})

test_that("fql coefficients are similar to NB GLM", {
  data(quine, package = "MASS")
  fit_fql <- fql(Days ~ Eth + Sex, data = quine)
  fit_nb <- MASS::glm.nb(Days ~ Eth + Sex, data = quine)

  # Coefficients should be in the same ballpark (within 50%)
  ratio <- abs(coef(fit_fql) / coef(fit_nb))
  expect_true(all(ratio > 0.5 & ratio < 2.0))
})

test_that("fql rejects invalid inputs", {
  data(quine, package = "MASS")

  expect_error(fql("not a formula", data = quine))
  expect_error(fql(Days ~ Eth, data = "not a data frame"))
  expect_error(fql(Days ~ Eth, data = quine, tol = -1))
  expect_error(fql(Days ~ Eth, data = quine, max.iter = 0))
})

test_that("fql with verbose mode runs without error", {
  data(quine, package = "MASS")
  expect_message(
    fql(Days ~ Eth, data = quine, verbose = TRUE),
    "Step 1"
  )
})
