test_that("source code includes bugfix guards", {
  fql_src <- paste(readLines(test_path("..", "..", "R", "fql.R"), warn = FALSE), collapse = "\n")
  fql_log_src <- paste(readLines(test_path("..", "..", "R", "fql_log.R"), warn = FALSE), collapse = "\n")

  expect_false(grepl("model\\.matrix\\(mt, mf, contrasts\\)", fql_src))
  expect_false(grepl("model\\.matrix\\(mt, mf, contrasts\\)", fql_log_src))

  expect_true(grepl("eta<-x%\\*%bold", fql_src))
  expect_true(grepl("eta<-x%\\*%bold", fql_log_src))

  compact_fql <- gsub("\\s+", "", fql_src)
  compact_fql_log <- gsub("\\s+", "", fql_log_src)
  expect_false(grepl("eta<-x%\\*%b}", compact_fql, perl = TRUE))
  expect_false(grepl("eta<-x%\\*%b}", compact_fql_log, perl = TRUE))

  expect_true(grepl("log\\(\\(y-mu\\)\\^2\\+eps\\)", compact_fql))
  expect_true(grepl("pmax\\(v,eps\\)", compact_fql))
  expect_true(grepl("pmax\\(v,eps\\)", compact_fql_log))
})

test_that("fql handles factor terms and returns positive variance estimates", {
  set.seed(123)
  n <- 80
  dat <- data.frame(
    y = rpois(n, lambda = 10),
    x = rnorm(n),
    g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
  )

  fit <- fql(y ~ x + g, data = dat, na.action = na.omit)
  expect_true(is.list(fit))
  expect_true(all(is.finite(fit$estimation$coefficients)))
  expect_true(all(fit$v > 0))
})
