# The estimates maximize the observed-data log-likelihood of Section 4.2 of the
# paper, written here independently with explicit loops over subjects and visits.
loglik_direct <- function(theta, d, K) {
  b1 <- theta[1:5]
  b2 <- theta[6:10]
  g1 <- theta[11:12]
  g2 <- theta[13:15]
  pr <- function(y, eta) {
    p <- 1 / (1 + exp(-eta))
    if (y == 1) p else 1 - p
  }
  ll <- 0
  for (i in unique(d$id)) {
    di <- d[d$id == i, ]
    Y1 <- Y2 <- rep(NA, K + 1)
    Y1[di$visit + 1] <- di$y1
    Y2[di$visit + 1] <- di$y2
    x <- di$x[1]
    a <- di$trt[1]
    if (!is.na(Y1[1])) ll <- ll + log(pr(Y1[1], g1[1] + g1[2] * x))
    if (!is.na(Y1[1]) && !is.na(Y2[1])) {
      ll <- ll + log(pr(Y2[1], g2[1] + g2[2] * Y1[1] + g2[3] * x))
    }
    for (k in seq_len(K)) {
      m1 <- !is.na(Y1[k + 1]) && !is.na(Y1[k]) && !is.na(Y2[k])       # m'_{i,1,k}
      m2 <- !is.na(Y2[k + 1]) && !is.na(Y2[k]) && !is.na(Y1[k + 1])   # m'_{i,2,k}
      if (m1) {
        ll <- ll + log(pr(Y1[k + 1], b1[1] + b1[2] * Y1[k] + b1[3] * Y2[k] + b1[4] * a +
                            b1[5] * x))
      }
      if (m2) {
        ll <- ll + log(pr(Y2[k + 1], b2[1] + b2[2] * Y1[k + 1] + b2[3] * Y2[k] + b2[4] * a +
                            b2[5] * x))
      }
    }
  }
  ll
}

test_that("the fit maximizes the missing-data likelihood of the paper", {
  d <- simulate_biscal(n = 150, K = 4, case = 2, missing = TRUE, seed = 21)
  set.seed(2)
  d$y2[d$visit == 0][sample(150, 8)] <- NA        # some missing baseline values
  d <- d[-sample(nrow(d), 40), ]                  # and some absent visits
  fit <- fit_binary(d)
  th <- coef(fit)
  expect_equal(fit$loglik, unname(loglik_direct(th, d, 4)), tolerance = 1e-10)
  h <- 1e-5
  score <- vapply(seq_along(th), function(i) {
    e <- replace(numeric(length(th)), i, h)
    (loglik_direct(th + e, d, 4) - loglik_direct(th - e, d, 4)) / (2 * h)
  }, numeric(1))
  expect_lt(max(abs(score)), 1e-4)
})
