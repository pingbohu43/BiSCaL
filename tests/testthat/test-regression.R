# Regression tests against reference values computed with the authors' original
# research code (see data-raw/regression-fixtures.R).
fx <- readRDS(test_path("fixtures", "original-code.rds"))

test_that("biscal() reproduces the authors' simulation code (complete data)", {
  s <- fx$sim_complete
  fit <- biscal(cbind(y1, y2) ~ x, data = s$long, id = "id", time = "visit",
                treatment = "trt")
  cf <- fit$coefficients
  # original order: gamma1, gamma2, beta1, beta2 = (1, Y2[k-1], Y1[k], trt, x)
  est <- c(cf$gamma1, cf$gamma2, cf$beta1, cf$beta2[c(1, 3, 2, 4, 5)])
  expect_equal(est, s$coef, tolerance = 1e-8)
  nm <- rownames(fit$vcov)
  ord <- c(grep("^gamma1", nm), grep("^gamma2", nm), grep("^beta1", nm),
           grep("^beta2", nm)[c(1, 3, 2, 4, 5)])
  expect_equal(unname(fit$vcov[ord, ord]), s$vcov, tolerance = 1e-8)

  # DEASTE averaged over the Monte Carlo draws (the original evaluated the
  # effects at K + 1, one step beyond the last visit)
  expect_warning(e1 <- easte(fit, newdata = data.frame(x = s$draws), K = s$K + 1),
                 "extrapolated")
  expect_equal(e1$test$estimate, s$deaste, tolerance = 1e-8)
  # its variance, with the derivative averaged over the subjects' covariates
  expect_warning(e2 <- easte(fit, newdata = data.frame(x = s$x), K = s$K + 1,
                             covariate_variance = FALSE), "extrapolated")
  expect_equal(e2$test$std.error^2, s$var, tolerance = 1e-6)
  # the 95% interval of the original code
  ci <- e1$test$estimate + c(-1, 1) * stats::qnorm(0.975) * e2$test$std.error
  expect_equal(ci, s$ci, tolerance = 1e-6)
})

test_that("true effects agree with the authors' simulation functions", {
  expect_equal(biscal_truth(K = 20, case = 1)[["DEASTE"]], fx$truth[["case1"]],
               tolerance = 1e-7)
  expect_equal(biscal_truth(K = 20, case = 3)[["DEASTE"]], fx$truth[["case3"]],
               tolerance = 1e-7)
})

test_that("effects agree with the authors' data-analysis functions (K = 64)", {
  r <- fx$realdata
  X <- rbind(r$xmin, r$xmax)
  v <- BiSCaL:::.ste_from_lgor(BiSCaL:::.lgor_arms(r$coefs, X, 64L, 0L, 0L))
  ste <- cbind(as.vector(v$STE1), as.vector(v$STE2))
  expect_equal(ste, unname(r$ste_paper_model), tolerance = 1e-8)
  # with these estimates P(Y_{2,0} = 1) is negligible, so the analysis code's
  # choice y_{2,0} = 0 in the model for Y_{1,1} does not change the effects
  expect_equal(ste, unname(r$ste_as_run), tolerance = 1e-8)
})
