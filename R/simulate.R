#' Coefficients of the simulation designs of the paper
#'
#' Returns the true coefficients of the three cases of the simulation studies
#' in Section 5.2 of the paper: binary responses, one standard normal covariate,
#' all coefficients equal to 0.5 except the treatment coefficients
#' \eqn{(\beta_{1,3,1}, \beta_{2,3,1})}, which are \eqn{(5, 0.1)} in case 1
#' (the treatment mainly affects the first response), \eqn{(0, 0)} in case 2
#' (no treatment effect) and \eqn{(0.1, 5)} in case 3 (the treatment mainly
#' affects the second response).
#'
#' @param case `1`, `2` or `3`.
#'
#' @return A list of one-column matrices `beta1`, `beta2`, `gamma1` and
#'   `gamma2`, laid out as the `coefficients` of a [biscal()] fit: terms in rows
#'   (intercept, the dummy-coded previous responses, the treatment `trt` and the
#'   covariate `x`, in the order of models (12)-(14) of the paper) and outcome
#'   categories in columns.
#'
#' @template ref
#' @seealso [simulate_biscal()], [biscal_truth()].
#'
#' @examples
#' biscal_coefs(1)
#'
#' @export
biscal_coefs <- function(case = 1) {
  case <- .check_case(case)
  trt <- switch(case, c(5, 0.1), c(0, 0), c(0.1, 5))
  mk <- function(v, terms) matrix(v, ncol = 1L, dimnames = list(terms, "1"))
  list(beta1 = mk(c(0.5, 0.5, 0.5, trt[1L], 0.5),
                  c("(Intercept)", "y1[k-1]=1", "y2[k-1]=1", "trt", "x")),
       beta2 = mk(c(0.5, 0.5, 0.5, trt[2L], 0.5),
                  c("(Intercept)", "y1[k]=1", "y2[k-1]=1", "trt", "x")),
       gamma1 = mk(c(0.5, 0.5), c("(Intercept)", "x")),
       gamma2 = mk(c(0.5, 0.5, 0.5), c("(Intercept)", "y1[0]=1", "x")))
}

.check_case <- function(case) {
  if (!(length(case) == 1L && !is.na(case) && as.character(case) %in% c("1", "2", "3"))) {
    .stopf("`case` must be 1, 2 or 3.")
  }
  as.integer(case)
}

# Check a list of coefficient matrices and return (L1, L2, p).
.check_coefs <- function(coefs) {
  nms <- c("beta1", "beta2", "gamma1", "gamma2")
  if (!is.list(coefs) || !all(nms %in% names(coefs))) {
    .stopf("`coefs` must be a list with matrices `beta1`, `beta2`, `gamma1` and `gamma2`.")
  }
  coefs <- lapply(coefs[nms], function(B) {
    B <- as.matrix(B)
    if (!is.numeric(B) || anyNA(B)) .stopf("The coefficients must be numeric and not missing.")
    B
  })
  L1 <- ncol(coefs$gamma1) + 1L
  L2 <- ncol(coefs$gamma2) + 1L
  p <- nrow(coefs$gamma1) - 1L
  ok <- p >= 0L && ncol(coefs$beta1) == L1 - 1L && ncol(coefs$beta2) == L2 - 1L &&
    nrow(coefs$gamma2) == L1 + p &&
    nrow(coefs$beta1) == L1 + L2 + p && nrow(coefs$beta2) == L1 + L2 + p
  if (!ok) {
    .stopf(paste0("The dimensions of the coefficient matrices in `coefs` are inconsistent; ",
                  "see biscal_coefs() for the layout."))
  }
  list(coefs = coefs, L1 = L1, L2 = L2, p = p)
}

#' Simulate longitudinal bivariate responses
#'
#' Generates data from the models of the paper: baseline responses from models
#' (13) and (14), then, at each visit \eqn{k = 1, \dots, K}, the first response
#' given the previous responses, the treatment and the covariates, and the
#' second response given the current first response, its own previous value,
#' the treatment and the covariates (model (12)). By default this is the design
#' of the simulation studies of the paper (Section 5.2): one standard normal
#' covariate, a treatment assigned with probability 0.5 independently of it,
#' and the coefficients of [biscal_coefs()].
#'
#' @section Missingness and censoring:
#' With `missing = TRUE` the data follow the design of Section 5.2.3 of the
#' paper. Subject \eqn{i} is followed up to visit
#' \eqn{C_i \sim \mathrm{Binomial}(K, \sigma(1.5 + x_i))} (no visit after
#' \eqn{C_i} is recorded), and at each visit \eqn{1 \le k \le C_i} the first and
#' second responses are missing independently with probabilities
#' \eqn{\sigma(x_i)} and \eqn{\sigma(x_i - 0.5)}, where \eqn{\sigma} is the
#' logistic function and \eqn{x_i} the first covariate. Baseline responses are
#' always observed.
#'
#' @param n Number of subjects.
#' @param K Number of follow-up visits; visits are numbered `0, 1, ..., K`.
#' @param case Simulation case of the paper (`1`, `2` or `3`), used when
#'   `coefs` is `NULL`.
#' @param missing Generate censoring and intermittent missingness as in the
#'   paper (see below)?
#' @param coefs Optional list of true coefficient matrices with the layout of
#'   [biscal_coefs()]; more columns give responses with more categories, and
#'   more covariate rows more covariates.
#' @param x Optional covariate values: a numeric vector or an `n`-row matrix
#'   with one column per covariate of `coefs`. By default independent standard
#'   normal variables.
#' @param seed Optional random seed. When given, the result is reproducible and
#'   the random number stream of the session is left unchanged.
#'
#' @return A data frame in long format, sorted by subject and visit, with
#'   columns `id`, `visit`, `trt` (0/1), the covariate(s) (`x` for a single
#'   covariate, otherwise `x1`, `x2`, ... or the column names of `x`), `y1`
#'   and `y2` (integer categories `0, 1, ...`; `NA` when missing). The
#'   attributes `coefs` and `censoring` hold the true coefficients and, with
#'   `missing = TRUE`, the last visit \eqn{C_i} of each subject.
#'
#' @template ref
#' @seealso [biscal_coefs()], [biscal_truth()], [biscal()].
#'
#' @examples
#' d <- simulate_biscal(n = 300, K = 5, case = 1, missing = TRUE, seed = 1)
#' head(d)
#' fit <- biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit",
#'               treatment = "trt")
#' coef(fit)
#'
#' @export
simulate_biscal <- function(n = 5000, K = 20, case = 1, missing = FALSE, coefs = NULL,
                            x = NULL, seed = NULL) {
  n <- .check_count(n, "n", min = 2L)
  K <- .check_count(K, "K", min = 1L)
  .check_flag(missing, "missing")
  seed <- .check_seed(seed)
  if (is.null(coefs)) coefs <- biscal_coefs(case)
  cc <- .check_coefs(coefs)
  p <- cc$p
  if (!is.null(seed)) {
    old <- .get_rng_state()
    on.exit(.restore_rng_state(old), add = TRUE)
    set.seed(seed)
  }
  if (is.null(x)) {
    X <- matrix(stats::rnorm(n * p), n, p)
    colnames(X) <- if (p == 1L) "x" else if (p > 1L) paste0("x", seq_len(p)) else NULL
  } else {
    X <- as.matrix(x)
    if (!is.numeric(X) || nrow(X) != n || ncol(X) != p || anyNA(X)) {
      .stopf("`x` must be numeric without missing values, with %d rows and %d column(s).", n, p)
    }
    if (is.null(colnames(X))) {
      colnames(X) <- if (p == 1L) "x" else paste0("x", seq_len(p))
    }
  }
  if (missing && p < 1L) .stopf("`missing = TRUE` requires at least one covariate.")
  A <- stats::rbinom(n, 1L, 0.5)
  Y <- .sim_responses(cc$coefs, X, A, K)

  out <- data.frame(id = rep(seq_len(n), each = K + 1L),
                    visit = rep(0:K, times = n),
                    trt = rep(A, each = K + 1L))
  for (j in seq_len(ncol(X))) out[[colnames(X)[j]]] <- rep(X[, j], each = K + 1L)
  out$y1 <- as.vector(t(Y$y1))
  out$y2 <- as.vector(t(Y$y2))
  cens <- NULL
  if (missing) {
    x1 <- X[, 1L]
    cens <- stats::rbinom(n, K, stats::plogis(1.5 + x1))
    m1 <- matrix(stats::runif(n * K) < stats::plogis(x1), n, K)
    m2 <- matrix(stats::runif(n * K) < stats::plogis(x1 - 0.5), n, K)
    y1 <- Y$y1
    y2 <- Y$y2
    y1[, -1L][m1] <- NA
    y2[, -1L][m2] <- NA
    out$y1 <- as.vector(t(y1))
    out$y2 <- as.vector(t(y2))
    out <- out[out$visit <= rep(cens, each = K + 1L), , drop = FALSE]
    rownames(out) <- NULL
  }
  attr(out, "coefs") <- cc$coefs
  attr(out, "censoring") <- cens
  out
}

# Responses from the models, one row per subject and one column per visit
# 0, ..., K (integer categories).
.sim_responses <- function(coefs, X, A, K) {
  n <- nrow(X)
  L1 <- ncol(coefs$gamma1) + 1L
  L2 <- ncol(coefs$gamma2) + 1L
  g1 <- coefs$gamma1
  g2 <- coefs$gamma2
  b1 <- coefs$beta1
  b2 <- coefs$beta2
  row_a <- L1 + L2
  lin <- function(B, D, xfirst) {
    eta <- .add_row(.xpart(B, X, xfirst), B[1L, ])
    if (!is.null(D) && ncol(D)) eta <- eta + D %*% B[1L + seq_len(ncol(D)), , drop = FALSE]
    eta
  }
  y1 <- matrix(NA_integer_, n, K + 1L)
  y2 <- matrix(NA_integer_, n, K + 1L)
  y1[, 1L] <- .rcat(.softmax0(lin(g1, NULL, 2L)))
  y2[, 1L] <- .rcat(.softmax0(lin(g2, .dummy(y1[, 1L], L1), L1 + 1L)))
  for (k in seq_len(K)) {
    D <- cbind(.dummy(y1[, k], L1), .dummy(y2[, k], L2), A)
    y1[, k + 1L] <- .rcat(.softmax0(lin(b1, D, row_a + 1L)))
    D <- cbind(.dummy(y1[, k + 1L], L1), .dummy(y2[, k], L2), A)
    y2[, k + 1L] <- .rcat(.softmax0(lin(b2, D, row_a + 1L)))
  }
  list(y1 = y1, y2 = y2)
}

#' True effects of a simulation design
#'
#' Computes the true separable and total treatment effects of a simulation
#' design with a single standard normal covariate (such as the designs of the
#' paper), either at given covariate values or averaged over the covariate
#' distribution, the targets of [aste()] and [easte()]. Population-level values
#' are computed by numerical integration.
#'
#' @param K End time point.
#' @param case,coefs The design, as in [simulate_biscal()].
#' @param y1,y2 Thresholds of the generalized odds ratio, as codes
#'   `0, 1, ...` of the response categories (not the highest).
#' @param x Optional numeric vector of covariate values. If `NULL` (default),
#'   population-level values under \eqn{X \sim N(0, 1)} are returned.
#'
#' @return With `x = NULL`, a named vector with the true `ESTE1`, `ESTE2`,
#'   `ETTE`, `EASTE1`, `EASTE2` and `DEASTE`. Otherwise a data frame with the
#'   true `STE1`, `STE2`, `TTE`, `ASTE1`, `ASTE2` and `DASTE` at each value of
#'   `x`.
#'
#' @template ref
#' @seealso [simulate_biscal()], [biscal_coefs()].
#'
#' @examples
#' # true DEASTE of the three cases (Table 5 of the paper)
#' sapply(1:3, function(case) biscal_truth(K = 20, case = case)[["DEASTE"]])
#' biscal_truth(K = 20, case = 1, x = c(-1, 0, 1))
#'
#' @export
biscal_truth <- function(K = 20, case = 1, coefs = NULL, y1 = 0, y2 = 0, x = NULL) {
  K <- .check_count(K, "K", min = 1L)
  if (is.null(coefs)) coefs <- biscal_coefs(case)
  cc <- .check_coefs(coefs)
  if (cc$p != 1L) .stopf("biscal_truth() supports designs with a single covariate.")
  c1 <- .check_count(y1, "y1", min = 0L)
  c2 <- .check_count(y2, "y2", min = 0L)
  if (c1 > cc$L1 - 2L || c2 > cc$L2 - 2L) {
    .stopf("`y1` and `y2` must be codes of categories below the highest one.")
  }
  eff <- function(xv) {
    v <- .ste_from_lgor(.lgor_arms(cc$coefs, matrix(xv, ncol = 1L), K, c1, c2))
    s1 <- as.vector(v$STE1)
    s2 <- as.vector(v$STE2)
    cbind(STE1 = s1, STE2 = s2, TTE = as.vector(v$TTE), ASTE1 = abs(s1),
          ASTE2 = abs(s2), DASTE = abs(s1) - abs(s2))
  }
  if (!is.null(x)) {
    if (!is.numeric(x) || anyNA(x)) .stopf("`x` must be a numeric vector.")
    return(data.frame(x = x, eff(x)))
  }
  vals <- vapply(1:6, function(e) {
    stats::integrate(function(t) eff(t)[, e] * stats::dnorm(t), -10, 10,
                     rel.tol = 1e-10, abs.tol = 1e-14, subdivisions = 1000L)$value
  }, numeric(1))
  stats::setNames(vals, .effect_names_pop)
}
