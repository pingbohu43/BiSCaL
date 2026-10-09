# Shared small data sets for the tests.

# The design of the paper with moderate treatment coefficients (in case 1 of
# the paper, the coefficient 5 separates the treated subjects in small samples).
coefs_moderate <- function() {
  co <- biscal_coefs(1)
  co$beta1["trt", 1] <- 1
  co$beta2["trt", 1] <- 0.3
  co
}

small_binary <- function(n = 300, K = 4, missing = FALSE, seed = 1) {
  simulate_biscal(n = n, K = K, coefs = coefs_moderate(), missing = missing, seed = seed)
}

# Run `expr`, returning its warnings (muffled).
collect_warnings <- function(expr) {
  w <- character(0)
  withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  w
}

coefs_3cat <- function() {
  mk <- function(B, terms) {
    matrix(B, ncol = 2L, dimnames = list(terms, c("1", "2")))
  }
  list(beta1 = mk(c(-0.5, 1.2, 2.0, 0.3, 0.5, 0.8, 0.4, -1.0, 0.8, 2.5, 0.2, 0.9, 0.2, 0.3),
                  c("(Intercept)", "y1[k-1]=1", "y1[k-1]=2", "y2[k-1]=1", "y2[k-1]=2", "trt", "x")),
       beta2 = mk(c(-0.3, 0.5, 0.9, 1.5, 2.2, -0.4, 0.2, -0.8, 0.7, 1.4, 1.0, 2.8, 0.6, -0.2),
                  c("(Intercept)", "y1[k]=1", "y1[k]=2", "y2[k-1]=1", "y2[k-1]=2", "trt", "x")),
       gamma1 = mk(c(0.2, 0.5, -0.3, 0.4), c("(Intercept)", "x")),
       gamma2 = mk(c(0.1, 0.8, 1.2, 0.3, -0.5, 0.4, 1.5, -0.2),
                   c("(Intercept)", "y1[0]=1", "y1[0]=2", "x")))
}

fit_binary <- function(d) {
  biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit", treatment = "trt")
}
