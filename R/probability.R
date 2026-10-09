# Counterfactual joint distribution of the two responses at the end time point,
# computed by the recursion of Theorem 2 / Theorem 3 of the paper:
#   P(Y_{1,k}(a1, a2) = y1, Y_{2,k}(a1, a2) = y2 | X = x)
#     = sum_{y1', y2'} P1(y1 | y1', y2', a1, x) P2(y2 | y2', y1, a2, x)
#                      P(Y_{1,k-1}(a1, a2) = y1', Y_{2,k-1}(a1, a2) = y2' | X = x),
# starting from P10(y1 | x) P20(y2 | y1, x) at k = 0 (equations (10)-(11)).
# The joint distribution is stored as a list of L1 * L2 vectors (one element
# per covariate profile); state (y1, y2) = (u, v) is element u + L1 * v + 1.

# Covariate contribution X %*% B[first:(first + p - 1), ] (an n x m matrix).
.xpart <- function(B, X, first) {
  p <- ncol(X)
  if (p == 0L) return(matrix(0, nrow(X), ncol(B)))
  X %*% B[first + seq_len(p) - 1L, , drop = FALSE]
}

# Add the row vector v to every row of the matrix M.
.add_row <- function(M, v) M + rep(v, each = nrow(M))

# Probabilities used by the recursion for the covariate profiles in X:
# the distribution at k = 0 (pi0, a list of S vectors) and, for a = 0, 1, the
# transition probabilities P1[[a + 1]][[s]][[u + 1]] = P(Y1 = u | previous
# state s = (y1', y2')) and P2[[a + 1]][[s]][[v + 1]] = P(Y2 = v | s = (current
# y1, previous y2')), each a vector over the profiles.
.cf_parts <- function(coefs, X) {
  n <- nrow(X)
  g1 <- coefs$gamma1
  g2 <- coefs$gamma2
  b1 <- coefs$beta1
  b2 <- coefs$beta2
  L1 <- ncol(g1) + 1L
  L2 <- ncol(g2) + 1L
  S <- L1 * L2
  su <- (seq_len(S) - 1L) %% L1
  sv <- (seq_len(S) - 1L) %/% L1

  P10 <- .softmax0(.add_row(.xpart(g1, X, 2L), g1[1L, ]))
  base20 <- .add_row(.xpart(g2, X, L1 + 1L), g2[1L, ])
  P20 <- lapply(seq_len(L1) - 1L, function(u) {
    .softmax0(if (u == 0L) base20 else .add_row(base20, g2[1L + u, ]))
  })
  pi0 <- lapply(seq_len(S), function(s) P10[, su[s] + 1L] * P20[[su[s] + 1L]][, sv[s] + 1L])

  row_a <- L1 + L2       # row of the treatment coefficient in beta1 and beta2
  base1 <- .add_row(.xpart(b1, X, row_a + 1L), b1[1L, ])
  base2 <- .add_row(.xpart(b2, X, row_a + 1L), b2[1L, ])
  offsets <- function(B, s, m) {
    off <- numeric(m)
    if (su[s] > 0L) off <- off + B[1L + su[s], ]
    if (sv[s] > 0L) off <- off + B[L1 + sv[s], ]
    off
  }
  cols <- function(P) lapply(seq_len(ncol(P)), function(j) P[, j])
  P1 <- P2 <- vector("list", 2L)
  for (a in 0:1) {
    lin1 <- if (a == 1L) .add_row(base1, b1[row_a, ]) else base1
    lin2 <- if (a == 1L) .add_row(base2, b2[row_a, ]) else base2
    P1[[a + 1L]] <- lapply(seq_len(S), function(s) {
      cols(.softmax0(.add_row(lin1, offsets(b1, s, L1 - 1L))))
    })
    P2[[a + 1L]] <- lapply(seq_len(S), function(s) {
      cols(.softmax0(.add_row(lin2, offsets(b2, s, L2 - 1L))))
    })
  }
  list(pi0 = pi0, P1 = P1, P2 = P2, L1 = L1, L2 = L2, S = S, su = su, sv = sv)
}

# One step of the recursion: first the first response moves from y1' to y1
# given (y1', y2'), then the second from y2' to y2 given (y1, y2').
.cf_step <- function(pi, P1a, P2a, L1, L2) {
  S <- L1 * L2
  q <- vector("list", S)
  for (v in seq_len(L2) - 1L) {
    prev <- v * L1 + seq_len(L1)          # states (y1' = 0, ..., L1 - 1; y2' = v)
    for (u in seq_len(L1)) {
      acc <- pi[[prev[1L]]] * P1a[[prev[1L]]][[u]]
      for (s0 in prev[-1L]) acc <- acc + pi[[s0]] * P1a[[s0]][[u]]
      q[[v * L1 + u]] <- acc
    }
  }
  out <- vector("list", S)
  for (u in seq_len(L1)) {
    prev <- u + L1 * (seq_len(L2) - 1L)   # states (y1 = u - 1; y2' = 0, ..., L2 - 1)
    for (v in seq_len(L2)) {
      acc <- q[[prev[1L]]] * P2a[[prev[1L]]][[v]]
      for (s0 in prev[-1L]) acc <- acc + q[[s0]] * P2a[[s0]][[v]]
      out[[u + L1 * (v - 1L)]] <- acc
    }
  }
  out
}

# Log generalized odds ratio (equation (2)) from a joint distribution pi, for
# thresholds c1 and c2 (codes of y1 and y2).
.lgor_from_pi <- function(pi, su, sv, c1, c2) {
  lo1 <- su <= c1
  lo2 <- sv <= c2
  pll <- Reduce(`+`, pi[lo1 & lo2])
  phh <- Reduce(`+`, pi[!lo1 & !lo2])
  plh <- Reduce(`+`, pi[lo1 & !lo2])
  phl <- Reduce(`+`, pi[!lo1 & lo2])
  log(pll) + log(phh) - log(plh) - log(phl)
}

# The four treatment arms in the order used throughout: (a1, a2) = (1, 1),
# (1, 0), (0, 1), (0, 0).
.arms <- rbind(c(1L, 1L), c(1L, 0L), c(0L, 1L), c(0L, 0L))

# log GOR_k(y1, y2, a1, a2, x) for every row of X, every k in Ks (sorted,
# unique) and the four arms: an n x length(Ks) x 4 array.
.lgor_arms <- function(coefs, X, Ks, c1, c2) {
  parts <- .cf_parts(coefs, X)
  n <- nrow(X)
  nK <- length(Ks)
  out <- array(NA_real_, c(n, nK, 4L))
  for (r in seq_len(4L)) {
    P1a <- parts$P1[[.arms[r, 1L] + 1L]]
    P2a <- parts$P2[[.arms[r, 2L] + 1L]]
    pi <- parts$pi0
    j <- 1L
    for (k in seq_len(Ks[nK])) {
      pi <- .cf_step(pi, P1a, P2a, parts$L1, parts$L2)
      if (k == Ks[j]) {
        out[, j, r] <- .lgor_from_pi(pi, parts$su, parts$sv, c1, c2)
        j <- j + 1L
      }
    }
  }
  out
}

# Joint distributions for one arm at the time points Ks: n x L1 x L2 x nK.
.cf_joint <- function(coefs, X, a1, a2, Ks) {
  parts <- .cf_parts(coefs, X)
  n <- nrow(X)
  nK <- length(Ks)
  out <- array(NA_real_, c(n, parts$L1, parts$L2, nK))
  pi <- parts$pi0
  j <- 1L
  if (Ks[1L] == 0L) {
    out[, , , 1L] <- unlist(pi, use.names = FALSE)
    j <- 2L
  }
  if (j <= nK) {
    for (k in seq_len(Ks[nK])) {
      pi <- .cf_step(pi, parts$P1[[a1 + 1L]], parts$P2[[a2 + 1L]], parts$L1, parts$L2)
      if (k == Ks[j]) {
        out[, , , j] <- unlist(pi, use.names = FALSE)
        j <- j + 1L
      }
    }
  }
  out
}

# Covariate matrix (without intercept) for a data frame of covariate profiles.
.newdata_X <- function(object, newdata, required = TRUE) {
  p <- ncol(object$X)
  if (is.null(newdata)) {
    if (p == 0L) return(matrix(numeric(0), nrow = 1L, ncol = 0L))
    if (!required) return(NULL)
    .stopf(paste0("`newdata` is required: a data frame of covariate values with one row ",
                  "per covariate profile and columns %s."),
           paste0("'", object$cov_vars, "'", collapse = ", "))
  }
  if (is.list(newdata) && !is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  if (!is.data.frame(newdata)) .stopf("`newdata` must be a data frame.")
  if (!nrow(newdata)) .stopf("`newdata` has no rows.")
  if (p == 0L) return(matrix(numeric(0), nrow = nrow(newdata), ncol = 0L))
  miss <- setdiff(object$cov_vars, names(newdata))
  if (length(miss)) {
    .stopf("`newdata` lacks the covariate(s) %s.", paste0("'", miss, "'", collapse = ", "))
  }
  mf <- stats::model.frame(object$terms, newdata[, object$cov_vars, drop = FALSE],
                           na.action = stats::na.pass, xlev = object$xlevels)
  if (!all(stats::complete.cases(mf))) .stopf("`newdata` contains missing covariate values.")
  X <- stats::model.matrix(object$terms, mf, contrasts.arg = object$contrasts)
  X <- X[, -1L, drop = FALSE]
  if (!identical(colnames(X), colnames(object$X))) {
    .stopf("The covariates in `newdata` do not match those of the fitted model.")
  }
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  rownames(X) <- NULL
  X
}

.resolve_K <- function(K, object, allow_zero = FALSE) {
  if (is.null(K)) return(object$K)
  lo <- if (allow_zero) 0 else 1
  if (!is.numeric(K) || !length(K) || anyNA(K) || any(!is.finite(K)) ||
      any(abs(K - round(K)) > 1e-8) || any(K < lo)) {
    .stopf("`K` must contain whole numbers >= %d (the end time point(s) of interest).", lo)
  }
  K <- sort(unique(as.integer(round(K))))
  if (any(K > object$K)) {
    .warnf(paste0("End time(s) %s lie beyond the last visit of the data (%d); the effects ",
                  "there are extrapolated with the fitted transition models."),
           paste(K[K > object$K], collapse = ", "), object$K)
  }
  K
}

.resolve_threshold <- function(y, lev, name, resp) {
  if (is.null(y)) return(0L)
  if (length(y) != 1L || is.na(y)) .stopf("`%s` must be a single category.", name)
  code <- match(as.character(y), lev) - 1L
  if (is.na(code)) {
    .stopf("`%s` must be one of the categories of response '%s': %s.", name, resp,
           paste(lev, collapse = ", "))
  }
  if (code == length(lev) - 1L) {
    .stopf("`%s` cannot be the highest category ('%s') of response '%s'.", name,
           lev[length(lev)], resp)
  }
  code
}

.check_biscal <- function(object) {
  if (!inherits(object, "biscal")) .stopf("`object` must be a fitted model returned by biscal().")
  invisible(object)
}

#' Counterfactual joint distribution of the two responses
#'
#' Computes \eqn{P\{Y_{1,k}(a_1, a_2) = y_1, Y_{2,k}(a_1, a_2) = y_2 \mid X = x\}},
#' the joint distribution of the two responses at visit \eqn{k} had the
#' treatment components been set to \eqn{(a_1, a_2)}, through the recursive
#' formula of Theorems 2 and 3 of the paper, from the fitted (or given) models.
#' The combinations \eqn{(1, 0)} and \eqn{(0, 1)} are counterfactual: they are
#' identified under the conditions of the paper.
#'
#' @param object A fitted model returned by [biscal()].
#' @param newdata A data frame of covariate values, one row per covariate
#'   profile \eqn{x}. May be omitted when the model has no covariates.
#' @param a1,a2 Values (0 or 1) of the treatment components acting on the
#'   first and the second response.
#' @param K Visit(s) \eqn{k} at which the distribution is computed; defaults to
#'   the last visit of the data. `0` gives the baseline distribution. Visits
#'   after the last visit of the data are extrapolated with the fitted
#'   transition models (with a warning).
#'
#' @return An array with dimensions (profile, category of the first response,
#'   category of the second response, visit), with dimnames.
#'
#' @template ref
#' @seealso [gor()] for the generalized odds ratios.
#'
#' @examples
#' fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
#'               time = "visit", treatment = "trt")
#' pr <- cf_prob(fit, newdata = data.frame(x1 = 0, x2 = 1), a1 = 1, a2 = 0)
#' pr[1, , , 1]          # joint distribution at the last visit
#' sum(pr[1, , , 1])     # sums to one
#'
#' @export
cf_prob <- function(object, newdata = NULL, a1, a2, K = NULL) {
  .check_biscal(object)
  for (nm in c("a1", "a2")) {
    v <- get(nm)
    if (!(length(v) == 1L && v %in% c(0, 1))) .stopf("`%s` must be 0 or 1.", nm)
  }
  X <- .newdata_X(object, newdata)
  Ks <- .resolve_K(K, object, allow_zero = TRUE)
  out <- .cf_joint(object$coefficients, X, as.integer(a1), as.integer(a2), Ks)
  dimnames(out) <- list(profile = seq_len(nrow(X)),
                        y1 = object$levels[[1L]], y2 = object$levels[[2L]],
                        K = Ks)
  names(dimnames(out))[2:3] <- object$response
  out
}

#' Generalized odds ratios of the two responses
#'
#' Computes the generalized odds ratio \eqn{GOR_k(y_1, y_2, a_1, a_2, x)} of
#' equation (2) of the paper: the odds ratio of the \eqn{2 \times 2} table
#' obtained by dichotomizing the counterfactual responses
#' \eqn{Y_{1,k}(a_1, a_2)} at \eqn{y_1} and \eqn{Y_{2,k}(a_1, a_2)} at
#' \eqn{y_2}, given \eqn{X = x},
#' \deqn{\frac{P(Y_1 \le y_1, Y_2 \le y_2) P(Y_1 > y_1, Y_2 > y_2)}{P(Y_1 \le y_1, Y_2 > y_2) P(Y_1 > y_1, Y_2 \le y_2)}.}{P(Y1 <= y1, Y2 <= y2) P(Y1 > y1, Y2 > y2) / {P(Y1 <= y1, Y2 > y2) P(Y1 > y1, Y2 <= y2)}.}
#' Values above one indicate that the two responses tend to fall on the same
#' side of their thresholds.
#'
#' @inheritParams cf_prob
#' @param y1,y2 Thresholds: categories of the first and second response (any
#'   category except the highest). Default to the lowest category, which is
#'   the only choice for binary responses.
#' @param K Visit(s) \eqn{k \ge 1} at which the ratios are computed; defaults
#'   to the last visit of the data. Visits after the last visit of the data are
#'   extrapolated with the fitted transition models (with a warning).
#' @param a1,a2 Values (0 or 1) of the treatment components; vectors of the
#'   same length give several arms at once.
#'
#' @return A data frame with columns `profile`, `K`, `a1`, `a2` and `gor`.
#'
#' @template ref
#' @seealso [aste()] for the separable treatment effects built from these
#'   ratios.
#'
#' @examples
#' fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
#'               time = "visit", treatment = "trt")
#' gor(fit, newdata = data.frame(x1 = 0, x2 = 1), a1 = c(1, 1, 0, 0),
#'     a2 = c(1, 0, 1, 0))
#'
#' @export
gor <- function(object, newdata = NULL, a1 = 1, a2 = 1, y1 = NULL, y2 = NULL, K = NULL) {
  .check_biscal(object)
  if (length(a1) != length(a2) || !length(a1) || !all(c(a1, a2) %in% c(0, 1))) {
    .stopf("`a1` and `a2` must be vectors of 0s and 1s of the same length.")
  }
  X <- .newdata_X(object, newdata)
  Ks <- .resolve_K(K, object)
  c1 <- .resolve_threshold(y1, object$levels[[1L]], "y1", object$response[1L])
  c2 <- .resolve_threshold(y2, object$levels[[2L]], "y2", object$response[2L])
  lg <- .lgor_arms(object$coefficients, X, Ks, c1, c2)
  n <- nrow(X)
  out <- lapply(seq_along(a1), function(i) {
    r <- which(.arms[, 1L] == a1[i] & .arms[, 2L] == a2[i])
    data.frame(profile = rep(seq_len(n), times = length(Ks)),
               K = rep(Ks, each = n), a1 = as.integer(a1[i]), a2 = as.integer(a2[i]),
               gor = exp(as.vector(lg[, , r])))
  })
  out <- do.call(rbind, out)
  out <- out[order(out$profile, out$K), , drop = FALSE]
  rownames(out) <- NULL
  out
}
