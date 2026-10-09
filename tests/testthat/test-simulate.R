test_that("biscal_coefs() gives the designs of the paper", {
  c1 <- biscal_coefs(1)
  expect_named(c1, c("beta1", "beta2", "gamma1", "gamma2"))
  expect_equal(c1$beta1["trt", 1], 5)
  expect_equal(c1$beta2["trt", 1], 0.1)
  expect_equal(biscal_coefs(2)$beta1["trt", 1], 0)
  expect_equal(biscal_coefs(3)$beta2["trt", 1], 5)
  expect_equal(unname(c1$gamma2[, 1]), rep(0.5, 3))
  expect_error(biscal_coefs(4), "1, 2 or 3")
})

test_that("simulate_biscal() returns long data with the requested design", {
  d <- simulate_biscal(n = 200, K = 6, case = 2, seed = 3)
  expect_named(d, c("id", "visit", "trt", "x", "y1", "y2"))
  expect_identical(nrow(d), 200L * 7L)
  expect_false(anyNA(d))
  expect_true(all(d$y1 %in% 0:1) && all(d$y2 %in% 0:1))
  expect_identical(attr(d, "coefs"), biscal_coefs(2))
  # covariate and treatment constant within subjects
  expect_true(all(tapply(d$x, d$id, function(v) length(unique(v))) == 1))
  expect_true(all(tapply(d$trt, d$id, function(v) length(unique(v))) == 1))
})

test_that("simulation is reproducible and leaves the RNG state unchanged", {
  set.seed(10)
  before <- .Random.seed
  a <- simulate_biscal(n = 50, K = 3, seed = 7)
  expect_identical(.Random.seed, before)
  b <- simulate_biscal(n = 50, K = 3, seed = 7)
  expect_identical(a, b)
  expect_false(identical(a, simulate_biscal(n = 50, K = 3, seed = 8)))
})

test_that("missingness and censoring follow Section 5.2.3 of the paper", {
  d <- simulate_biscal(n = 4000, K = 20, case = 1, missing = TRUE, seed = 11)
  cens <- attr(d, "censoring")
  expect_length(cens, 4000L)
  last <- tapply(d$visit, d$id, max)
  expect_identical(as.integer(last), as.integer(cens))
  # mean censoring time and mean missingness rates reported in the paper
  # (missingness rates averaged over subjects followed beyond baseline)
  expect_equal(mean(cens), 15.571, tolerance = 0.02)
  fu <- d[d$visit > 0, ]
  expect_equal(mean(tapply(is.na(fu$y1), fu$id, mean)), 0.5, tolerance = 0.03)
  expect_equal(mean(tapply(is.na(fu$y2), fu$id, mean)), 0.398, tolerance = 0.03)
  expect_false(anyNA(d$y1[d$visit == 0]))
})

test_that("simulate_biscal() accepts custom coefficients and covariates", {
  d <- simulate_biscal(n = 300, K = 2, coefs = coefs_3cat(), seed = 1)
  expect_true(all(d$y1 %in% 0:2) && all(d$y2 %in% 0:2))
  co <- biscal_coefs(1)
  co$gamma1 <- rbind(co$gamma1, x2 = 0.2)
  co$gamma2 <- rbind(co$gamma2, x2 = 0.2)
  co$beta1 <- rbind(co$beta1, x2 = 0.2)
  co$beta2 <- rbind(co$beta2, x2 = 0.2)
  X <- cbind(age = rnorm(100), sex = rbinom(100, 1, 0.5))
  d2 <- simulate_biscal(n = 100, K = 2, coefs = co, x = X, seed = 1)
  expect_named(d2, c("id", "visit", "trt", "age", "sex", "y1", "y2"))
  expect_error(simulate_biscal(n = 100, K = 2, coefs = co, x = X[, 1]), "column")
  co$beta1 <- co$beta1[-1, , drop = FALSE]
  expect_error(simulate_biscal(n = 100, K = 2, coefs = co), "inconsistent")
  expect_error(simulate_biscal(n = 1), "whole number >= 2")
})

test_that("estimates from simulated data are close to the truth", {
  d <- simulate_biscal(n = 3000, K = 5, coefs = coefs_moderate(), seed = 2)
  fit <- fit_binary(d)
  z <- (coef(fit) - BiSCaL:::.flatten(coefs_moderate())) / sqrt(diag(vcov(fit)))
  expect_lt(max(abs(z)), 4)
})

test_that("biscal_truth() gives the true effects", {
  t2 <- biscal_truth(K = 5, case = 2)
  expect_equal(unname(t2), rep(0, 6), tolerance = 1e-10)
  t1 <- biscal_truth(K = 5, case = 1)
  expect_named(t1, c("ESTE1", "ESTE2", "ETTE", "EASTE1", "EASTE2", "DEASTE"))
  expect_equal(t1[["DEASTE"]], t1[["EASTE1"]] - t1[["EASTE2"]], tolerance = 1e-8)
  expect_equal(t1[["ESTE1"]] + t1[["ESTE2"]], 2 * t1[["ETTE"]], tolerance = 1e-8)
  ind <- biscal_truth(K = 5, case = 1, x = c(-1, 0, 1))
  expect_identical(nrow(ind), 3L)
  expect_equal(ind$DASTE, ind$ASTE1 - ind$ASTE2)
  # Monte Carlo check of the population value
  set.seed(1)
  mc <- mean(biscal_truth(K = 5, case = 1, x = rnorm(4000))$DASTE)
  expect_equal(t1[["DEASTE"]], mc, tolerance = 0.05)
  expect_error(biscal_truth(K = 5, coefs = coefs_3cat(), y1 = 2), "below the highest")
})
