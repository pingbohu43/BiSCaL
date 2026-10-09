test_that("cf_prob() returns probability distributions", {
  fit <- fit_binary(small_binary(missing = TRUE))
  nd <- data.frame(x = c(-1, 0, 2))
  pr <- cf_prob(fit, newdata = nd, a1 = 1, a2 = 0, K = c(0, 2, 4))
  expect_identical(dim(pr), c(3L, 2L, 2L, 3L))
  expect_identical(names(dimnames(pr)), c("profile", "y1", "y2", "K"))
  expect_equal(apply(pr, c(1, 4), sum), matrix(1, 3, 3), tolerance = 1e-12, ignore_attr = TRUE)
  # at baseline the treatment components play no role
  pr0 <- cf_prob(fit, newdata = nd, a1 = 0, a2 = 1, K = 0)
  expect_equal(pr[, , , 1], pr0[, , , 1], tolerance = 1e-14)
  expect_error(cf_prob(fit, newdata = nd, a1 = 2, a2 = 0), "0 or 1")
  expect_error(cf_prob(fit, a1 = 1, a2 = 1), "newdata")
})

test_that("the recursion agrees with a Monte Carlo simulation of the chain", {
  co <- coefs_3cat()
  x0 <- 0.7
  K <- 3
  set.seed(5)
  N <- 2e5
  X <- matrix(x0, N, 1)
  arm <- c(1, 0)
  y1 <- BiSCaL:::.rcat(BiSCaL:::.softmax0(BiSCaL:::.add_row(X %*% co$gamma1[2, , drop = FALSE],
                                                          co$gamma1[1, ])))
  eta2 <- BiSCaL:::.add_row(X %*% co$gamma2[4, , drop = FALSE], co$gamma2[1, ]) +
    BiSCaL:::.dummy(y1, 3) %*% co$gamma2[2:3, ]
  y2 <- BiSCaL:::.rcat(BiSCaL:::.softmax0(eta2))
  for (k in seq_len(K)) {
    e1 <- BiSCaL:::.add_row(X %*% co$beta1[7, , drop = FALSE], co$beta1[1, ]) +
      cbind(BiSCaL:::.dummy(y1, 3), BiSCaL:::.dummy(y2, 3), arm[1]) %*% co$beta1[2:6, ]
    y1n <- BiSCaL:::.rcat(BiSCaL:::.softmax0(e1))
    e2 <- BiSCaL:::.add_row(X %*% co$beta2[7, , drop = FALSE], co$beta2[1, ]) +
      cbind(BiSCaL:::.dummy(y1n, 3), BiSCaL:::.dummy(y2, 3), arm[2]) %*% co$beta2[2:6, ]
    y2 <- BiSCaL:::.rcat(BiSCaL:::.softmax0(e2))
    y1 <- y1n
  }
  mc <- unclass(table(factor(y1, 0:2), factor(y2, 0:2))) / N
  rec <- BiSCaL:::.cf_joint(co, matrix(x0, 1), arm[1], arm[2], K)[1, , , 1]
  expect_lt(max(abs(rec - mc)), 4 * sqrt(0.25 / N))
})

test_that("gor() matches the joint distribution", {
  fit <- fit_binary(small_binary())
  nd <- data.frame(x = 0.4)
  pr <- cf_prob(fit, newdata = nd, a1 = 1, a2 = 1, K = 3)[1, , , 1]
  g <- gor(fit, newdata = nd, a1 = 1, a2 = 1, K = 3)
  expect_equal(g$gor, pr[1, 1] * pr[2, 2] / (pr[1, 2] * pr[2, 1]), tolerance = 1e-12)
  g4 <- gor(fit, newdata = data.frame(x = c(0, 1)), a1 = c(1, 1, 0, 0),
            a2 = c(1, 0, 1, 0), K = c(2, 4))
  expect_identical(nrow(g4), 16L)
  expect_true(all(g4$gor > 0))
  expect_error(gor(fit, newdata = nd, a1 = c(1, 0), a2 = 1), "same length")
})

test_that("thresholds collapse ordinal responses as in equation (2)", {
  co <- coefs_3cat()
  d <- simulate_biscal(n = 1500, K = 3, coefs = co, seed = 8)
  fit <- biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit", treatment = "trt")
  nd <- data.frame(x = -0.2)
  pr <- cf_prob(fit, newdata = nd, a1 = 0, a2 = 1)[1, , , 1]
  for (c1 in 0:1) for (c2 in 0:1) {
    lo1 <- seq_len(3) <= c1 + 1
    lo2 <- seq_len(3) <= c2 + 1
    want <- sum(pr[lo1, lo2]) * sum(pr[!lo1, !lo2]) /
      (sum(pr[lo1, !lo2]) * sum(pr[!lo1, lo2]))
    got <- gor(fit, newdata = nd, a1 = 0, a2 = 1, y1 = c1, y2 = c2)$gor
    expect_equal(got, want, tolerance = 1e-12)
  }
  expect_error(gor(fit, newdata = nd, y1 = 2), "highest")
  expect_error(gor(fit, newdata = nd, y2 = "z"), "one of the categories")
})
