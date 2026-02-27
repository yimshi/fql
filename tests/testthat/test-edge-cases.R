test_that("fql handles many covariates", {
  data(quine, package = "MASS")
  fit <- fql(Days ~ Eth + Sex + Age + Lrn + Eth:Sex, data = quine)

  expect_s3_class(fit, "fql")
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("fql handles data with zeros", {
  set.seed(123)
  n <- 100
  x <- rnorm(n)
  # Generate data with some zeros
  mu <- exp(1 + 0.5 * x)
  y <- rpois(n, lambda = mu)
  y[sample(n, 20)] <- 0  # add extra zeros
  dat <- data.frame(y = y, x = x)

  fit <- fql(y ~ x, data = dat)
  expect_s3_class(fit, "fql")
})

test_that("fql_multi basic functionality", {
  set.seed(42)
  n <- 80
  otu <- matrix(MASS::rnegbin(n * 3, mu = 10, theta = 2), nrow = n)
  colnames(otu) <- paste0("OTU", 1:3)
  dat <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)))

  results <- fql_multi(otu, ~ group, data = dat,
                       prevalence.filter = 0, verbose = FALSE)

  expect_true(is.data.frame(results))
  expect_true("taxon" %in% names(results))
  expect_true("p_adjusted" %in% names(results))
  expect_true(all(results$p_adjusted >= results$p_value, na.rm = TRUE))
})

test_that("fql_multi prevalence filter works", {
  set.seed(42)
  n <- 80
  otu <- matrix(MASS::rnegbin(n * 3, mu = 10, theta = 2), nrow = n)
  # Make one OTU very sparse
  otu[, 3] <- 0
  otu[1:5, 3] <- c(1, 2, 1, 3, 1)
  colnames(otu) <- paste0("OTU", 1:3)
  dat <- data.frame(group = factor(rep(c("A", "B"), each = n / 2)))

  results <- fql_multi(otu, ~ group, data = dat,
                       prevalence.filter = 0.25, verbose = FALSE)

  # OTU3 should be filtered out (prevalence < 25%)
  expect_true(!"OTU3" %in% results$taxon)
})

test_that("fql_multi rejects invalid inputs", {
  expect_error(fql_multi("not a matrix", ~ x, data = data.frame(x = 1)))
})

test_that("sim_microbiome dataset loads correctly", {
  data(sim_microbiome, package = "fql")

  expect_true(is.list(sim_microbiome))
  expect_true("otu_table" %in% names(sim_microbiome))
  expect_true("sample_data" %in% names(sim_microbiome))
  expect_true("taxa_info" %in% names(sim_microbiome))

  expect_equal(nrow(sim_microbiome$otu_table), 200)
  expect_equal(ncol(sim_microbiome$otu_table), 10)
  expect_equal(nrow(sim_microbiome$sample_data), 200)
})
