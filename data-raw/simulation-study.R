# Reproduces the simulation studies of Section 5.2 of the paper (Tables 4-7)
# with BiSCaL and saves the results used by vignette("simulation"):
# inst/extdata/simulation-study.rds. Run from the package root, e.g.
#   Rscript data-raw/simulation-study.R
# (about 2 seconds per replicate and core; 600 replicates in total).
devtools::load_all(quiet = TRUE)
n <- 5000
K <- 20
n_draws <- 10000
R <- 100
cores <- max(1L, parallel::detectCores())

one_rep <- function(i, case, missing) {
  seed <- 100000L * missing + 1000L * case + i
  d <- simulate_biscal(n = n, K = K, case = case, missing = missing, seed = seed)
  warn <- character(0)
  fit <- withCallingHandlers(
    biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit", treatment = "trt"),
    warning = function(w) {
      warn <<- c(warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  # the expectation over X ~ N(0, 1) by Monte Carlo draws, as in the paper
  set.seed(seed + 500000L)
  draws <- data.frame(x = rnorm(n_draws))
  e <- easte(fit, newdata = draws, K = K)
  list(coef = coef(fit), se = sqrt(diag(vcov(fit))), deaste = e$test$estimate,
       se_deaste = e$test$std.error, p = e$test$p.value, warnings = length(warn))
}

scenarios <- expand.grid(case = 1:3, missing = c(FALSE, TRUE))
results <- lapply(seq_len(nrow(scenarios)), function(s) {
  case <- scenarios$case[s]
  missing <- scenarios$missing[s]
  t0 <- proc.time()[["elapsed"]]
  reps <- parallel::mclapply(seq_len(R), one_rep, case = case, missing = missing,
                             mc.cores = cores)
  message(sprintf("case %d, missing = %s: %.0f s", case, missing,
                  proc.time()[["elapsed"]] - t0))
  list(case = case, missing = missing,
       coef = do.call(rbind, lapply(reps, `[[`, "coef")),
       se = do.call(rbind, lapply(reps, `[[`, "se")),
       deaste = vapply(reps, `[[`, numeric(1), "deaste"),
       se_deaste = vapply(reps, `[[`, numeric(1), "se_deaste"),
       p = vapply(reps, `[[`, numeric(1), "p"),
       warnings = vapply(reps, `[[`, numeric(1), "warnings"),
       true_coef = .flatten(biscal_coefs(case)),
       true_deaste = biscal_truth(K = K, case = case)[["DEASTE"]])
})
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(settings = list(n = n, K = K, n_draws = n_draws, R = R), results = results),
        "inst/extdata/simulation-study.rds", version = 2, compress = "xz")
