# Maximum likelihood for one baseline-category logistic regression model,
#   log P(Y = d | z) / P(Y = 0 | z) = z' b_d,   d = 1, ..., L - 1,
# the form of models (12)-(14) of the paper. Binary outcomes (L = 2) are
# fitted by stats::glm.fit(), exactly as in the authors' code; outcomes with
# more categories by Newton-Raphson with step halving.
#
# Returns the q x (L - 1) coefficient matrix, its covariance matrix (inverse
# observed information, coefficients stacked by outcome category), the
# log-likelihood, the number of observations and convergence information.

.biscal_control <- function(control) {
  if (!is.list(control)) .stopf("`control` must be a list.")
  ok <- c("epsilon", "maxit")
  bad <- setdiff(names(control), ok)
  if (length(bad)) .stopf("Unknown element(s) of `control`: %s.", paste(bad, collapse = ", "))
  out <- list(epsilon = 1e-8, maxit = 25L)
  out[names(control)] <- control
  .check_number(out$epsilon, "control$epsilon", positive = TRUE)
  out$maxit <- .check_count(out$maxit, "control$maxit")
  out
}

.fit_mlogit <- function(Z, y, L, lev, label, control) {
  n <- nrow(Z)
  q <- ncol(Z)
  m <- L - 1L
  if (n == 0L) {
    .stopf(paste0("No observations are available to fit the %s.%s"), label,
           if (grepl("^transition", label)) paste0(
             " A transition needs the responses at two consecutive visits k - 1 and k; ",
             "check that the visits are numbered 0, 1, 2, ... without gaps.") else "")
  }
  tab <- tabulate(y + 1L, nbins = L)
  if (any(tab == 0L)) {
    .stopf("Category '%s' never occurs as the outcome of the %s, which cannot be fitted.",
           lev[which(tab == 0L)[1L]], label)
  }
  qz <- qr(Z)
  if (qz$rank < q) {
    alias <- colnames(Z)[qz$pivot[(qz$rank + 1L):q]]
    .stopf(paste0("The %s cannot be fitted: the data carry no information to separate ",
                  "the effect of %s from the other terms (for example, a category or a ",
                  "treatment group that is absent among the contributing observations)."),
           label, paste0("'", alias, "'", collapse = ", "))
  }
  res <- if (m == 1L) {
    .fit_binary(Z, y, label, control)
  } else {
    .fit_multinomial(Z, y, L, tab, label, control)
  }
  se <- sqrt(pmax(diag(res$vcov), 0))
  big <- which(se > 100)
  if (length(big)) {
    terms <- if (m == 1L) colnames(Z) else
      paste0(rep(colnames(Z), m), " (category ", rep(lev[-1L], each = q), ")")
    .warnf(paste0("In the %s: very large standard errors (%s) indicate (quasi-)complete ",
                  "separation; the corresponding fitted probabilities are practically 0 or 1."),
           label, paste0("'", terms[big], "'", collapse = ", "))
  }
  res
}

.fit_binary <- function(Z, y, label, control) {
  q <- ncol(Z)
  msgs <- character(0)
  fit <- withCallingHandlers(
    stats::glm.fit(Z, y, family = stats::binomial(),
                   control = stats::glm.control(epsilon = control$epsilon,
                                                maxit = control$maxit)),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (length(msgs)) .warnf("In the %s: %s", label, paste(unique(msgs), collapse = "; "))
  Qr <- fit$qr
  p1 <- seq_len(fit$rank)
  covu <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  V <- matrix(NA_real_, q, q)
  V[Qr$pivot[p1], Qr$pivot[p1]] <- covu
  B <- matrix(unname(fit$coefficients), ncol = 1L)
  list(coef = B, vcov = V, loglik = -fit$deviance / 2, nobs = length(y),
       converged = isTRUE(fit$converged), iter = fit$iter)
}

# Log-likelihood of a baseline-category logit model (log-sum-exp form).
.mlogit_loglik <- function(Z, y, B) {
  eta <- cbind(0, Z %*% B)
  mx <- apply(eta, 1L, max)
  lse <- mx + log(rowSums(exp(eta - mx)))
  sum(eta[cbind(seq_len(nrow(Z)), y + 1L)] - lse)
}

# Observed (= expected) information matrix, coefficients stacked by category.
.mlogit_info <- function(Z, P) {
  q <- ncol(Z)
  m <- ncol(P)
  H <- matrix(0, q * m, q * m)
  for (d in seq_len(m)) {
    id <- (d - 1L) * q + seq_len(q)
    for (e in d:m) {
      ie <- (e - 1L) * q + seq_len(q)
      w <- if (d == e) P[, d] * (1 - P[, d]) else -P[, d] * P[, e]
      blk <- crossprod(Z, Z * w)
      H[id, ie] <- blk
      if (e != d) H[ie, id] <- t(blk)
    }
  }
  H
}

.fit_multinomial <- function(Z, y, L, tab, label, control) {
  n <- nrow(Z)
  q <- ncol(Z)
  m <- L - 1L
  Y <- .dummy(y, L)
  B <- matrix(0, q, m)
  B[1L, ] <- log(tab[-1L] / tab[1L])
  ll <- .mlogit_loglik(Z, y, B)
  maxit <- max(control$maxit, 100L)
  converged <- FALSE
  iter <- 0L
  for (iter in seq_len(maxit)) {
    P <- .softmax0(Z %*% B)[, -1L, drop = FALSE]
    g <- as.vector(crossprod(Z, Y - P))
    H <- .mlogit_info(Z, P)
    step <- tryCatch(solve(H, g), error = function(e) NULL)
    if (is.null(step)) break
    t <- 1
    repeat {
      Bn <- B + t * matrix(step, q, m)
      lln <- .mlogit_loglik(Z, y, Bn)
      if (is.finite(lln) && lln >= ll - 1e-12 * abs(ll)) break
      t <- t / 2
      if (t < 1e-10) break
    }
    if (!is.finite(lln) || t < 1e-10) break
    change <- abs(lln - ll) / (abs(lln) + 0.1)
    B <- Bn
    ll <- lln
    if (change < control$epsilon) {
      converged <- TRUE
      break
    }
  }
  if (!converged) .warnf("In the %s: the algorithm did not converge.", label)
  Pfull <- .softmax0(Z %*% B)
  if (any(Pfull < 10 * .Machine$double.eps)) {
    .warnf("In the %s: fitted probabilities numerically 0 occurred.", label)
  }
  H <- .mlogit_info(Z, Pfull[, -1L, drop = FALSE])
  V <- tryCatch(chol2inv(chol(H)), error = function(e) {
    tryCatch(solve(H), error = function(e2) matrix(NA_real_, q * m, q * m))
  })
  if (anyNA(V)) .warnf("In the %s: the information matrix is singular.", label)
  list(coef = B, vcov = V, loglik = ll, nobs = n, converged = converged, iter = iter)
}
