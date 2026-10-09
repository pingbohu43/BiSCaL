test_that("biscal() fits the four models and returns a consistent object", {
  d <- small_binary(missing = TRUE)
  fit <- fit_binary(d)
  expect_s3_class(fit, "biscal")
  expect_named(fit$coefficients, c("beta1", "beta2", "gamma1", "gamma2"))
  expect_identical(rownames(fit$coefficients$beta1),
                   c("(Intercept)", "y1[k-1]=1", "y2[k-1]=1", "trt", "x"))
  expect_identical(rownames(fit$coefficients$beta2),
                   c("(Intercept)", "y1[k]=1", "y2[k-1]=1", "trt", "x"))
  expect_identical(rownames(fit$coefficients$gamma2), c("(Intercept)", "y1[0]=1", "x"))
  cf <- coef(fit)
  expect_length(cf, 15L)
  expect_identical(names(cf), rownames(vcov(fit)))
  expect_identical(names(cf)[1:2], c("beta1:(Intercept)", "beta1:y1[k-1]=1"))
  expect_true(isSymmetric(unname(vcov(fit))))
  expect_true(all(diag(vcov(fit)) > 0))
  # block-diagonal: the four likelihood factors share no parameters
  expect_equal(unname(vcov(fit)[1:5, 6:15]), matrix(0, 5, 10))
  expect_identical(nobs(fit), 300L)
  expect_identical(fit$K, 4L)
  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_identical(attr(ll, "df"), 15L)
  expect_equal(as.numeric(ll), sum(fit$loglik_models))
  expect_true(all(fit$converged))
})

test_that("binary models agree with glm()", {
  d <- small_binary()
  fit <- fit_binary(d)
  w <- reshape(d[, c("id", "visit", "y1", "y2")], idvar = "id", timevar = "visit",
               direction = "wide")
  base <- unique(d[, c("id", "trt", "x")])
  g1 <- glm(w$y1.0 ~ base$x, family = binomial)
  expect_equal(unname(fit$coefficients$gamma1[, 1]), unname(coef(g1)), tolerance = 1e-10)
  idx <- grep("^gamma1", rownames(vcov(fit)))
  expect_equal(unname(vcov(fit)[idx, idx]), unname(vcov(g1)), tolerance = 1e-10)
  # transition model for y1 by stacking the visits
  tr <- do.call(rbind, lapply(1:4, function(k) {
    data.frame(y = w[[paste0("y1.", k)]], l1 = w[[paste0("y1.", k - 1)]],
               l2 = w[[paste0("y2.", k - 1)]], trt = base$trt, x = base$x)
  }))
  g3 <- glm(y ~ l1 + l2 + trt + x, family = binomial, data = tr)
  expect_equal(unname(fit$coefficients$beta1[, 1]), unname(coef(g3)), tolerance = 1e-10)
  expect_equal(as.numeric(fit$loglik_models[["beta1"]]), as.numeric(logLik(g3)),
               tolerance = 1e-10)
})

test_that("the likelihood uses the transitions required by m'_{i,j,k}", {
  d <- small_binary(n = 200, K = 3, missing = TRUE, seed = 4)
  d$y2[d$id == 1 & d$visit == 0] <- NA
  d <- d[!(d$id == 2 & d$visit == 1), ]
  fit <- fit_binary(d)
  Y1 <- Y2 <- matrix(NA, 200, 4)
  Y1[cbind(d$id, d$visit + 1)] <- d$y1
  Y2[cbind(d$id, d$visit + 1)] <- d$y2
  n1 <- sum(sapply(1:3, function(k) sum(!is.na(Y1[, k + 1]) & !is.na(Y1[, k]) & !is.na(Y2[, k]))))
  n2 <- sum(sapply(1:3, function(k) sum(!is.na(Y2[, k + 1]) & !is.na(Y2[, k]) & !is.na(Y1[, k + 1]))))
  expect_equal(unname(fit$nobs_models), c(n1, n2, sum(!is.na(Y1[, 1])),
                                          sum(!is.na(Y1[, 1]) & !is.na(Y2[, 1]))))
})

test_that("rows can be in any order and identifiers of any type", {
  d <- small_binary(missing = TRUE)
  fit <- fit_binary(d)
  set.seed(9)
  d2 <- d[sample(nrow(d)), ]
  d2$id <- paste0("S", d2$id)
  fit2 <- fit_binary(d2)
  expect_equal(coef(fit2), coef(fit), tolerance = 1e-10)
})

test_that("responses and treatment may be factors or logical", {
  d <- small_binary()
  fit <- fit_binary(d)
  d2 <- d
  d2$y1 <- factor(ifelse(d$y1 == 1, "high", "low"), levels = c("low", "high"))
  d2$y2 <- d$y2 == 1
  d2$trt <- factor(ifelse(d$trt == 1, "drug", "placebo"), levels = c("placebo", "drug"))
  fit2 <- biscal(cbind(y1, y2) ~ x, data = d2, id = "id", time = "visit", treatment = "trt")
  expect_equal(unname(coef(fit2)), unname(coef(fit)), tolerance = 1e-10)
  expect_identical(fit2$levels, list(c("low", "high"), c("FALSE", "TRUE")))
  expect_identical(fit2$treatment_levels, c("placebo", "drug"))
  expect_identical(rownames(fit2$coefficients$beta1)[2], "y1[k-1]=high")
})

test_that("factor covariates, transformations and models without covariates work", {
  d <- small_binary(n = 400)
  d$g <- factor(c("a", "b", "c"))[(d$id %% 3) + 1]
  fit <- biscal(cbind(y1, y2) ~ x + I(x^2) + g, data = d, id = "id", time = "visit",
                treatment = "trt")
  expect_identical(rownames(fit$coefficients$gamma1),
                   c("(Intercept)", "x", "I(x^2)", "gb", "gc"))
  res <- aste(fit, newdata = data.frame(x = 0.3, g = "c"))
  expect_true(all(is.finite(res$effects$estimate)))
  expect_error(aste(fit, newdata = data.frame(x = 0.3, g = "z")), "new level")
  fit0 <- biscal(cbind(y1, y2) ~ 1, data = d, id = "id", time = "visit", treatment = "trt")
  expect_identical(rownames(fit0$coefficients$beta1),
                   c("(Intercept)", "y1[k-1]=1", "y2[k-1]=1", "trt"))
  expect_s3_class(aste(fit0), "biscal_aste")
  expect_s3_class(easte(fit0), "biscal_easte")
})

test_that("invalid input gives informative errors", {
  d <- small_binary()
  f <- cbind(y1, y2) ~ x
  expect_error(biscal(f, data = as.list(d), id = "id", time = "visit", treatment = "trt"),
               "data frame")
  expect_error(biscal(y1 ~ x, data = d, id = "id", time = "visit", treatment = "trt"),
               "cbind")
  expect_error(biscal(f, data = d, id = "ID", time = "visit", treatment = "trt"),
               "not in `data`")
  expect_error(biscal(cbind(y1, y2) ~ ., data = d, id = "id", time = "visit",
                      treatment = "trt"), "not supported")
  expect_error(biscal(cbind(y1, y2) ~ x - 1, data = d, id = "id", time = "visit",
                      treatment = "trt"), "intercept")
  expect_error(biscal(cbind(y1, y2) ~ x + trt, data = d, id = "id", time = "visit",
                      treatment = "trt"), "cannot be used as a covariate")
  d2 <- d
  d2$visit <- d2$visit + 1
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "baseline")
  d2 <- d
  d2$visit[2] <- 0.5
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "visit index")
  d2 <- rbind(d, d[1, ])
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "more than one row")
  d2 <- d
  d2$x[2] <- d2$x[2] + 1
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "varies within subjects")
  d2 <- d
  d2$trt[2] <- 1 - d2$trt[2]
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "constant within subjects")
  d2 <- d
  d2$trt <- d2$trt + 1
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "coded 0")
  d2 <- d
  d2$trt <- 0
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "Both treatment groups")
  d2 <- d
  d2$x[d2$id == 3] <- NA
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "missing for 1 subject")
  d2 <- d
  d2$y1 <- as.character(d2$y1)
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "ordered")
  d2 <- d
  d2$y2 <- 1
  expect_error(biscal(f, data = d2, id = "id", time = "visit", treatment = "trt"),
               "at least two")
  d2 <- d
  d2$x2 <- 2 * d2$x
  expect_error(biscal(cbind(y1, y2) ~ x + x2, data = d2, id = "id", time = "visit",
                      treatment = "trt"), "collinear")
  expect_error(biscal(f, data = d, id = "id", time = "visit", treatment = "trt",
                      control = list(tol = 1)), "Unknown")
})

test_that("unused factor levels are dropped with a warning", {
  d <- small_binary()
  d$y1 <- factor(d$y1, levels = c(0, 1, 2))
  expect_warning(fit <- fit_binary(d), "Unused categories")
  expect_identical(fit$levels[[1]], c("0", "1"))
})

test_that("categories absent from a model are reported", {
  d <- small_binary()
  # the first response is 0 at every follow-up visit: category 1 never occurs
  # as the outcome of its transition model
  d1 <- d
  d1$y1[d1$visit > 0] <- 0
  expect_error(fit_binary(d1), "never occurs")
  # category 1 occurs only at the last visit, so never as a previous value
  d2 <- d
  d2$y1[d2$visit < max(d2$visit)] <- 0
  expect_error(fit_binary(d2), "no information to separate the effect of 'y1\\[k-1\\]=1'")
})

test_that("printing and summaries work", {
  fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id", time = "visit",
                treatment = "trt")
  expect_output(print(fit), "Transition model for y1")
  s <- summary(fit)
  expect_s3_class(s, "summary.biscal")
  expect_named(s$coefficients, c("beta1", "beta2", "gamma1", "gamma2"))
  expect_identical(colnames(s$coefficients$beta1),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  expect_output(print(s), "Baseline model for y2")
  ci <- confint(fit, parm = "beta1:trt", level = 0.9)
  expect_equal(unname(ci[1, ]), unname(coef(fit)["beta1:trt"] +
                                         c(-1, 1) * qnorm(0.95) * sqrt(vcov(fit)["beta1:trt", "beta1:trt"])))
  expect_identical(nrow(confint(fit)), length(coef(fit)))
  expect_error(confint(fit, parm = "nope"), "Unknown")
})

test_that("data-dependent transformations use the parameters of the fitted data", {
  fit <- biscal(cbind(y1, y2) ~ poly(x1, 2) + x2, data = simtrial,
                id = "id", time = "visit", treatment = "trt")
  nd <- data.frame(x1 = c(-1, 0.5, 2), x2 = c(0, 1, 1))
  all3 <- aste(fit, newdata = nd)$effects
  one <- aste(fit, newdata = nd[2, ])$effects
  expect_equal(one$estimate, all3$estimate[all3$profile == 2], tolerance = 1e-12)
  # same fit with the covariate standardized by hand
  base <- unique(simtrial[, c("id", "x1")])
  m <- mean(base$x1)
  sdv <- sd(base$x1)
  d2 <- simtrial
  d2$z <- (d2$x1 - m) / sdv
  fit2 <- biscal(cbind(y1, y2) ~ z + x2, data = d2, id = "id", time = "visit",
                 treatment = "trt")
  fit1 <- biscal(cbind(y1, y2) ~ scale(x1) + x2, data = simtrial, id = "id",
                 time = "visit", treatment = "trt")
  expect_equal(unname(coef(fit1)), unname(coef(fit2)), tolerance = 1e-10)
  e1 <- aste(fit1, newdata = data.frame(x1 = 1.3, x2 = 0))$effects$estimate
  e2 <- aste(fit2, newdata = data.frame(z = (1.3 - m) / sdv, x2 = 0))$effects$estimate
  expect_equal(e1, e2, tolerance = 1e-10)
})

test_that("unused covariate levels are dropped and offsets rejected", {
  d <- simtrial
  d$g <- factor(ifelse(d$x2 == 1, "b", "a"), levels = c("a", "b", "c"))
  fit <- biscal(cbind(y1, y2) ~ x1 + g, data = d, id = "id", time = "visit",
                treatment = "trt")
  expect_identical(fit$xlevels$g, c("a", "b"))
  expect_error(aste(fit, newdata = data.frame(x1 = 0, g = "c")), "new level")
  expect_error(biscal(cbind(y1, y2) ~ x1 + offset(x2), data = simtrial, id = "id",
                      time = "visit", treatment = "trt"), "offset")
})

test_that("visits that are not numbered consecutively are reported", {
  d <- small_binary()
  d$visit <- 2 * d$visit
  expect_error(fit_binary(d), "consecutive visits")
})

test_that("effects beyond the last visit come with a warning", {
  fit <- fit_binary(small_binary())
  expect_warning(aste(fit, newdata = data.frame(x = 0), K = 6), "extrapolated")
  expect_warning(gor(fit, newdata = data.frame(x = 0), K = c(2, 9)), "End time\\(s\\) 9")
  expect_no_warning(aste(fit, newdata = data.frame(x = 0), K = c(4, 4, 2)))
})

test_that("print() shows the treatment coding", {
  d <- small_binary()
  d$trt <- factor(ifelse(d$trt == 1, "drug", "placebo"), levels = c("placebo", "drug"))
  fit <- fit_binary(d)
  expect_output(print(fit), "control = placebo, treated = drug")
})
