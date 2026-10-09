test_that("aste() returns effects, tests and GORs consistent with Theorem 1", {
  fit <- fit_binary(small_binary(n = 500, K = 5, missing = TRUE))
  nd <- data.frame(x = c(-1, 0.5))
  res <- aste(fit, newdata = nd, K = c(3, 5))
  expect_s3_class(res, "biscal_aste")
  expect_identical(nrow(res$effects), 2L * 2L * 6L)
  expect_identical(levels(res$effects$effect), c("STE1", "STE2", "TTE", "ASTE1", "ASTE2", "DASTE"))
  w <- reshape(res$effects[, c("profile", "K", "effect", "estimate")],
               idvar = c("profile", "K"), timevar = "effect", direction = "wide")
  expect_equal(w$estimate.STE1 + w$estimate.STE2, 2 * w$estimate.TTE, tolerance = 1e-12)
  expect_equal(w$estimate.DASTE, abs(w$estimate.STE1) - abs(w$estimate.STE2))
  expect_true(all(res$effects$std.error > 0))
  # effects from the GORs
  g <- res$gor[res$gor$profile == 1 & res$gor$K == 5, ]
  lg <- log(g$gor)   # arms (1,1), (1,0), (0,1), (0,0)
  ste1 <- res$effects$estimate[res$effects$profile == 1 & res$effects$K == 5 &
                                 res$effects$effect == "STE1"]
  expect_equal(ste1, (lg[1] - lg[3]) + (lg[2] - lg[4]), tolerance = 1e-12)
  # a vector of end times gives the same values as separate calls
  r5 <- aste(fit, newdata = nd, K = 5)
  expect_equal(r5$effects[, c("estimate", "std.error")],
               res$effects[res$effects$K == 5, c("estimate", "std.error")],
               ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("delta-method standard errors match a direct computation", {
  fit <- fit_binary(small_binary(n = 400, K = 4))
  nd <- data.frame(x = 0.3)
  res <- aste(fit, newdata = nd)
  f <- function(theta) {
    co <- BiSCaL:::.unflatten(theta, fit$coefficients)
    v <- BiSCaL:::.ste_from_lgor(BiSCaL:::.lgor_arms(co, matrix(0.3, 1), 4L, 0L, 0L))
    c(v$STE1, v$STE2)
  }
  th <- coef(fit)
  h <- 1e-5
  J <- sapply(seq_along(th), function(i) {
    e <- replace(numeric(length(th)), i, h)
    (f(th + e) - f(th - e)) / (2 * h)
  })
  se <- sqrt(diag(J %*% vcov(fit) %*% t(J)))
  got <- res$effects$std.error[res$effects$effect %in% c("STE1", "STE2")]
  expect_equal(got, se, tolerance = 1e-5)
  s <- sign(f(th))
  gd <- s[1] * J[1, ] - s[2] * J[2, ]
  expect_equal(res$test$std.error, sqrt(drop(gd %*% vcov(fit) %*% gd)), tolerance = 1e-5)
})

test_that("p-values and intervals follow the requested alternative and level", {
  fit <- fit_binary(small_binary(n = 400, K = 4))
  nd <- data.frame(x = 0)
  lo <- aste(fit, newdata = nd, alternative = "less")$test
  gr <- aste(fit, newdata = nd, alternative = "greater")$test
  ts <- aste(fit, newdata = nd, alternative = "two.sided")$test
  expect_equal(lo$p.value, pnorm(lo$statistic))
  expect_equal(lo$p.value + gr$p.value, 1)
  expect_equal(ts$p.value, 2 * min(lo$p.value, gr$p.value))
  r90 <- aste(fit, newdata = nd, conf.level = 0.9)$effects
  expect_equal(r90$conf.high - r90$estimate, qnorm(0.95) * r90$std.error)
  expect_error(aste(fit, newdata = nd, conf.level = 1.2), "conf.level")
  expect_error(aste(fit, newdata = nd, K = 0), "whole numbers >= 1")
  expect_error(aste(fit), "newdata")
  expect_error(aste(fit, newdata = data.frame(z = 1)), "lacks the covariate")
  expect_error(aste(fit, newdata = data.frame(x = NA_real_)), "missing covariate")
  expect_error(aste(list()), "biscal")
})

test_that("easte() averages over covariates and adds their sampling variance", {
  fit <- fit_binary(small_binary(n = 400, K = 4, missing = TRUE))
  e <- easte(fit)
  expect_s3_class(e, "biscal_easte")
  expect_identical(e$n_rows, 400L)
  X <- unique(small_binary(n = 400, K = 4, missing = TRUE)[, c("id", "x")])
  ind <- aste(fit, newdata = data.frame(x = X$x))
  pop <- c(STE1 = "ESTE1", STE2 = "ESTE2", TTE = "ETTE", ASTE1 = "EASTE1",
           ASTE2 = "EASTE2", DASTE = "DEASTE")
  for (nm in names(pop)) {
    expect_equal(e$effects$estimate[e$effects$effect == pop[[nm]]],
                 mean(ind$effects$estimate[ind$effects$effect == nm]), tolerance = 1e-12)
  }
  h <- ind$effects$estimate[ind$effects$effect == "DASTE"]
  expect_equal(e$effects$mc_se[e$effects$effect == "DEASTE"], sd(h) / sqrt(length(h)),
               tolerance = 1e-10)
  e0 <- easte(fit, covariate_variance = FALSE)
  expect_equal(e0$effects$estimate, e$effects$estimate)
  expect_true(all(e$effects$std.error >= e0$effects$std.error))
  expect_equal(e$effects$std.error^2, e0$effects$std.error^2 + e$effects$mc_se^2,
               tolerance = 1e-10)
  # supplying the subjects' covariates as newdata excludes the covariate variance
  e2 <- easte(fit, newdata = data.frame(x = X$x))
  expect_equal(e2$effects$std.error, e0$effects$std.error, tolerance = 1e-12)
  expect_false(e2$covariate_variance)
  # computation in chunks of rows gives the same results
  op <- options(BiSCaL.chunk_budget = 2000)
  e3 <- easte(fit)
  a3 <- aste(fit, newdata = data.frame(x = X$x[1:70]))
  options(op)
  expect_equal(e3$effects[, -1], e$effects[, -1], tolerance = 1e-10)
  a1 <- aste(fit, newdata = data.frame(x = X$x[1:70]))
  expect_equal(a3$effects, a1$effects, tolerance = 1e-12)
  ek <- easte(fit, K = c(2, 4))
  expect_equal(ek$effects$estimate[ek$effects$K == 4], e$effects$estimate, tolerance = 1e-12)
})

test_that("effects vanish when the treatment has no effect", {
  co <- biscal_coefs(2)
  v <- BiSCaL:::.ste_from_lgor(BiSCaL:::.lgor_arms(co, matrix(c(-1, 0, 1)), 5L, 0L, 0L))
  expect_equal(as.vector(v$STE1), rep(0, 3), tolerance = 1e-12)
  expect_equal(as.vector(v$TTE), rep(0, 3), tolerance = 1e-12)
})

test_that("separable effects of ordinal responses can be estimated", {
  d <- simulate_biscal(n = 1500, K = 4, coefs = coefs_3cat(), missing = TRUE, seed = 2)
  fit <- biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit", treatment = "trt")
  expect_identical(dim(fit$coefficients$beta1), c(7L, 2L))
  expect_identical(names(coef(fit))[1], "beta1[1]:(Intercept)")
  res <- aste(fit, newdata = data.frame(x = 0), y1 = 1, y2 = 0)
  expect_identical(res$y1, "1")
  w <- res$effects$estimate
  expect_equal(w[1] + w[2], 2 * w[3], tolerance = 1e-12)
  expect_true(all(is.finite(res$effects$std.error)))
  expect_output(print(res), "y1 <= 1, y2 <= 0")
})

test_that("printing and conversion of effects work", {
  fit <- fit_binary(small_binary(n = 300, K = 3))
  res <- aste(fit, newdata = data.frame(x = c(-1, 1)))
  expect_output(print(res), "Profile 2 \\(x = 1\\)")
  expect_output(print(res), "H0: DASTE >= 0 vs H1: DASTE < 0")
  many <- aste(fit, newdata = data.frame(x = c(-1, 1)), K = 1:3)
  expect_output(print(many, max_print = 2), "Estimates of STE1, STE2, TTE and DASTE")
  expect_identical(as.data.frame(res), res$effects)
  expect_identical(as.data.frame(res, what = "test"), res$test)
  expect_identical(as.data.frame(res, what = "gor"), res$gor)
  e <- easte(fit, alternative = "two.sided")
  expect_output(print(e), "H0: DEASTE = 0 vs H1: DEASTE != 0")
  expect_identical(as.data.frame(e, what = "test"), e$test)
})
