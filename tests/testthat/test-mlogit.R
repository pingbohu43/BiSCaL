test_that("multinomial models agree with nnet::multinom()", {
  skip_if_not_installed("nnet")
  d <- simulate_biscal(n = 1500, K = 3, coefs = coefs_3cat(), missing = TRUE, seed = 3)
  fit <- biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit", treatment = "trt")
  des <- BiSCaL:::.build_designs(BiSCaL:::.prepare_data(cbind(y1, y2) ~ x, d, "id",
                                                        "visit", "trt"))
  for (nm in names(des)) {
    Z <- des[[nm]]$Z
    y <- factor(des[[nm]]$y)
    mm <- nnet::multinom(y ~ Z - 1, trace = FALSE, Hess = TRUE, reltol = 1e-14,
                         maxit = 1000)
    expect_equal(as.vector(fit$coefficients[[nm]]), as.vector(t(coef(mm))),
                 tolerance = 1e-5)
    idx <- grep(paste0("^", nm), rownames(fit$vcov))
    expect_equal(sqrt(diag(fit$vcov[idx, idx])), sqrt(diag(solve(mm$Hessian))),
                 tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fit$loglik_models[[nm]], as.numeric(logLik(mm)), tolerance = 1e-8)
  }
})

test_that("Newton-Raphson for two categories agrees with glm.fit()", {
  set.seed(4)
  Z <- cbind(1, matrix(rnorm(600), 300))
  y <- rbinom(300, 1, plogis(Z %*% c(0.3, 1, -0.5)))
  ctrl <- BiSCaL:::.biscal_control(list())
  nr <- BiSCaL:::.fit_multinomial(Z, y, 2L, tabulate(y + 1L, 2L), "test", ctrl)
  gl <- BiSCaL:::.fit_binary(Z, y, "test", ctrl)
  expect_equal(nr$coef, gl$coef, tolerance = 1e-8)
  # glm.fit() computes the covariance from the weights of its last iteration
  expect_equal(nr$vcov, gl$vcov, tolerance = 1e-4)
  expect_equal(nr$loglik, gl$loglik, tolerance = 1e-10)
})

test_that("separation is reported", {
  set.seed(5)
  ctrl <- BiSCaL:::.biscal_control(list())
  # complete separation: glm.fit() warns
  x <- rnorm(100)
  Z <- cbind("(Intercept)" = 1, x = x)
  w <- collect_warnings(BiSCaL:::.fit_mlogit(Z, as.integer(x > 0), 2L, c("0", "1"),
                                             "toy model", ctrl))
  expect_true(any(grepl("^In the toy model: glm.fit", w)))
  # quasi-complete separation: huge standard errors
  Z <- cbind("(Intercept)" = 1, g = c(rep(0, 50), rep(1, 50)))
  y <- c(rbinom(50, 1, 0.5), rep(1, 50))
  expect_warning(BiSCaL:::.fit_mlogit(Z, y, 2L, c("0", "1"), "toy model", ctrl),
                 "very large standard errors.*'g'")
})
