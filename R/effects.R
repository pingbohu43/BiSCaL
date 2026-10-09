# Separable and total treatment effects, delta-method standard errors and
# Wald tests at the individual level (given X = x, Section 3.2.1 of the paper)
# and at the population level (averaged over a covariate distribution,
# Section 3.2.2).

.effect_names_ind <- c("STE1", "STE2", "TTE", "ASTE1", "ASTE2", "DASTE")
.effect_names_pop <- c("ESTE1", "ESTE2", "ETTE", "EASTE1", "EASTE2", "DEASTE")

# Step of the central-difference derivatives: the default of pracma::grad(),
# which the authors' code used.
.deriv_step <- .Machine$double.eps^(1 / 3)

# STE1, STE2 and TTE (equations (3) and (4)) from the log GORs of the four arms,
# stored in the order (1, 1), (1, 0), (0, 1), (0, 0) along the third dimension
# of an n x nK x 4 array (values) or an n x nK x 4 x P array (derivatives).
.ste_from_lgor <- function(l) {
  if (length(dim(l)) == 3L) {
    pick <- function(r) l[, , r, drop = FALSE]
  } else {
    pick <- function(r) l[, , r, , drop = FALSE]
  }
  l11 <- pick(1L)
  l10 <- pick(2L)
  l01 <- pick(3L)
  l00 <- pick(4L)
  list(STE1 = (l11 - l01) + (l10 - l00),
       STE2 = (l11 - l10) + (l01 - l00),
       TTE = l11 - l00)
}

# Rows of X processed together, so that the arrays of derivatives stay small
# (at most about `budget` numbers; the option is used by the tests).
.chunk_size <- function(nK, P, budget = getOption("BiSCaL.chunk_budget", 4e6)) {
  max(1L, as.integer(budget %/% (4 * nK * (P + 1))))
}

# Values and central-difference derivatives (with respect to the stacked
# coefficients) of the log GORs of the four arms for the rows of X.
# Returns value (n x nK x 4) and jac (n x nK x 4 x P).
.lgor_with_jacobian <- function(coefs, X, Ks, c1, c2, step = .deriv_step) {
  theta <- .flatten(coefs)
  P <- length(theta)
  value <- .lgor_arms(coefs, X, Ks, c1, c2)
  jac <- array(0, c(dim(value), P))
  for (i in seq_len(P)) {
    tp <- theta
    tm <- theta
    tp[i] <- tp[i] + step
    tm[i] <- tm[i] - step
    jac[, , , i] <- (.lgor_arms(.unflatten(tp, coefs), X, Ks, c1, c2) -
                       .lgor_arms(.unflatten(tm, coefs), X, Ks, c1, c2)) / (2 * step)
  }
  list(value = value, jac = jac)
}

.check_finite_effects <- function(v) {
  if (any(!is.finite(v))) {
    .stopf(paste0("The effects could not be computed for some covariate values: the ",
                  "fitted probabilities of some cells of the 2 x 2 table are numerically ",
                  "zero. Use less extreme covariate values."))
  }
}

# Individual level: per profile and end time, the estimates of
# STE1, STE2, TTE, ASTE1, ASTE2, DASTE and their delta-method standard errors.
.effects_individual <- function(object, X, Ks, c1, c2) {
  coefs <- object$coefficients
  V <- object$vcov
  n <- nrow(X)
  nK <- length(Ks)
  P <- nrow(V)
  est <- se <- array(NA_real_, c(n, nK, 6L))
  lgor <- array(NA_real_, c(n, nK, 4L))
  cs <- .chunk_size(nK, P)
  for (first in seq.int(1L, n, by = cs)) {
    rows <- first:min(n, first + cs - 1L)
    J <- .lgor_with_jacobian(coefs, X[rows, , drop = FALSE], Ks, c1, c2)
    .check_finite_effects(J$value)
    lgor[rows, , ] <- J$value
    v <- .ste_from_lgor(J$value)
    g <- .ste_from_lgor(J$jac)
    s1 <- sign(v$STE1)
    s2 <- sign(v$STE2)
    vals <- list(v$STE1, v$STE2, v$TTE, abs(v$STE1), abs(v$STE2), abs(v$STE1) - abs(v$STE2))
    grads <- list(g$STE1, g$STE2, g$TTE,
                  sweep(g$STE1, 1:3, s1, `*`), sweep(g$STE2, 1:3, s2, `*`),
                  sweep(g$STE1, 1:3, s1, `*`) - sweep(g$STE2, 1:3, s2, `*`))
    for (e in 1:6) {
      est[rows, , e] <- vals[[e]][, , 1L]
      G <- matrix(grads[[e]], ncol = P)   # (rows x nK) x P
      se[rows, , e] <- sqrt(pmax(rowSums((G %*% V) * G), 0))
    }
  }
  list(est = est, se = se, lgor = lgor)
}

# Population level: averages over the rows of X of STE1, STE2, TTE, |STE1|,
# |STE2| and |STE1| - |STE2|, the derivatives of these averages, and the
# sample variance of the averaged quantities over the rows.
.effects_population <- function(object, X, Ks, c1, c2) {
  coefs <- object$coefficients
  n <- nrow(X)
  nK <- length(Ks)
  P <- nrow(object$vcov)
  sums <- matrix(0, nK, 6L)
  grad <- array(0, c(nK, 6L, P))
  cmean <- matrix(0, nK, 6L)
  cm2 <- matrix(0, nK, 6L)
  done <- 0L
  cs <- .chunk_size(nK, P)
  for (first in seq.int(1L, n, by = cs)) {
    rows <- first:min(n, first + cs - 1L)
    nb <- length(rows)
    J <- .lgor_with_jacobian(coefs, X[rows, , drop = FALSE], Ks, c1, c2)
    .check_finite_effects(J$value)
    v <- .ste_from_lgor(J$value)
    g <- .ste_from_lgor(J$jac)
    s1 <- v$STE1[, , 1L]
    s2 <- v$STE2[, , 1L]
    h <- list(s1, s2, v$TTE[, , 1L], abs(s1), abs(s2), abs(s1) - abs(s2))
    w1 <- sign(v$STE1)
    w2 <- sign(v$STE2)
    dg <- list(g$STE1, g$STE2, g$TTE, sweep(g$STE1, 1:3, w1, `*`),
               sweep(g$STE2, 1:3, w2, `*`),
               sweep(g$STE1, 1:3, w1, `*`) - sweep(g$STE2, 1:3, w2, `*`))
    for (e in 1:6) {
      he <- matrix(h[[e]], nb, nK)
      sums[, e] <- sums[, e] + colSums(he)
      grad[, e, ] <- grad[, e, ] + matrix(colSums(matrix(dg[[e]], nb)), nK, P)
      # pooled mean and sum of squared deviations (Chan et al.)
      bm <- colMeans(he)
      bm2 <- colSums((he - rep(bm, each = nb))^2)
      delta <- bm - cmean[, e]
      tot <- done + nb
      cm2[, e] <- cm2[, e] + bm2 + delta^2 * done * nb / tot
      cmean[, e] <- cmean[, e] + delta * nb / tot
    }
    done <- done + nb
  }
  list(est = sums / n, grad = grad / n,
       var_rows = if (n > 1L) cm2 / (n - 1L) else matrix(0, nK, 6L))
}

.p_value <- function(z, alternative) {
  switch(alternative,
         less = stats::pnorm(z),
         greater = stats::pnorm(z, lower.tail = FALSE),
         two.sided = 2 * stats::pnorm(-abs(z)))
}

.profile_labels <- function(object, newdata, n) {
  if (!length(object$cov_vars) || is.null(newdata)) return(rep("(no covariates)", n))
  parts <- lapply(object$cov_vars, function(v) {
    x <- newdata[[v]]
    val <- if (is.numeric(x)) format(signif(x, 4), trim = TRUE) else as.character(x)
    paste0(v, " = ", val)
  })
  do.call(paste, c(parts, list(sep = ", ")))
}

#' Separable treatment effects and tests given covariate values
#'
#' Estimates, for given covariate values \eqn{x}, the separable treatment
#' effects \eqn{STE_{1,K}} and \eqn{STE_{2,K}} on the two responses, the total
#' treatment effect \eqn{TTE_K}, the absolute separable treatment effects
#' \eqn{ASTE_{j,K} = |STE_{j,K}|} and their difference
#' \eqn{DASTE_K = ASTE_{1,K} - ASTE_{2,K}}, with delta-method standard errors,
#' and tests whether the treatment affects the second response more than the
#' first (Section 3.2.1 of the paper).
#'
#' @section Definitions:
#' With \eqn{GOR_K(a_1, a_2)} the generalized odds ratio of the responses at
#' the end time point \eqn{K} under treatment components \eqn{(a_1, a_2)} (see
#' [gor()]),
#' \deqn{TTE_K = \log\frac{GOR_K(1,1)}{GOR_K(0,0)},}{TTE_K = log{GOR_K(1,1) / GOR_K(0,0)},}
#' \deqn{STE_{1,K} = \log\frac{GOR_K(1,1)}{GOR_K(0,1)} + \log\frac{GOR_K(1,0)}{GOR_K(0,0)},}{STE_1K = log{GOR_K(1,1) / GOR_K(0,1)} + log{GOR_K(1,0) / GOR_K(0,0)},}
#' \deqn{STE_{2,K} = \log\frac{GOR_K(1,1)}{GOR_K(1,0)} + \log\frac{GOR_K(0,1)}{GOR_K(0,0)},}{STE_2K = log{GOR_K(1,1) / GOR_K(1,0)} + log{GOR_K(0,1) / GOR_K(0,0)},}
#' so that \eqn{STE_{1,K} + STE_{2,K} = 2\,TTE_K} (Theorem 1 of the paper).
#' \eqn{STE_{1,K}} changes the treatment component acting on the first
#' response with the other component held fixed, and \eqn{STE_{2,K}}
#' conversely.
#'
#' @section Inference:
#' Standard errors use the delta method with the covariance matrix of the
#' maximum likelihood estimators (equation (24) of the paper); the derivatives
#' are central differences with step \eqn{\epsilon^{1/3}} (\eqn{\epsilon} the
#' machine precision), as in the authors' implementation, and the absolute
#' values are differentiated through \eqn{\mathrm{sign}(STE_{j,K})}, a
#' subgradient. The default test is that of the paper,
#' \eqn{H_0: DASTE_K \ge 0} against \eqn{H_1: DASTE_K < 0} (the absolute
#' separable effect on the second response is larger), with p-value
#' \eqn{\Phi(z)}, \eqn{z = \widehat{DASTE}_K / \widehat{se}}. Confidence
#' intervals are Wald intervals.
#'
#' @param object A fitted model returned by [biscal()].
#' @param newdata A data frame of covariate values, one row per covariate
#'   profile \eqn{x} (for example, a row of the data, or the smallest or largest
#'   value of each covariate as in Table 2 of the paper). May be omitted when
#'   the model has no covariates.
#' @param y1,y2 Thresholds of the generalized odds ratio: categories of the
#'   first and second response (any category except the highest). Default to
#'   the lowest category, the only choice for binary responses.
#' @param K End time point(s): the visit(s) at which the effects are defined.
#'   Defaults to the last visit of the data; a vector gives effects at several
#'   visits (see [plot.biscal_aste()]). Visits after the last visit of the data
#'   are extrapolated with the fitted transition models (with a warning).
#' @param alternative Alternative hypothesis for the test of \eqn{DASTE_K}
#'   (\eqn{DEASTE_K} for [easte()]): `"less"` (default, as in the paper),
#'   `"greater"` or `"two.sided"`.
#' @param conf.level Confidence level of the Wald intervals.
#'
#' @return An object of class `"biscal_aste"`, a list with components
#'   \describe{
#'     \item{`effects`}{Data frame with one row per profile, end time and
#'       effect (`STE1`, `STE2`, `TTE`, `ASTE1`, `ASTE2`, `DASTE`): `estimate`,
#'       `std.error`, `conf.low`, `conf.high`.}
#'     \item{`test`}{Data frame with the test of \eqn{DASTE_K} for each profile
#'       and end time: `estimate`, `std.error`, `statistic` (\eqn{z}) and
#'       `p.value`.}
#'     \item{`gor`}{Data frame of the estimated generalized odds ratios of the
#'       four arms.}
#'   }
#'   and the settings used.
#'
#' @template ref
#' @seealso [easte()] for the population-level effects; [gor()];
#'   [plot.biscal_aste()].
#'
#' @examples
#' fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
#'               time = "visit", treatment = "trt")
#' # effects for two covariate profiles
#' profiles <- data.frame(x1 = c(-1, 1), x2 = c(0, 1))
#' res <- aste(fit, newdata = profiles)
#' res
#' as.data.frame(res)
#'
#' @export
aste <- function(object, newdata = NULL, y1 = NULL, y2 = NULL, K = NULL,
                 alternative = c("less", "greater", "two.sided"), conf.level = 0.95) {
  .check_biscal(object)
  alternative <- match.arg(alternative)
  .check_conf_level(conf.level)
  X <- .newdata_X(object, newdata)
  Ks <- .resolve_K(K, object)
  c1 <- .resolve_threshold(y1, object$levels[[1L]], "y1", object$response[1L])
  c2 <- .resolve_threshold(y2, object$levels[[2L]], "y2", object$response[2L])
  res <- .effects_individual(object, X, Ks, c1, c2)
  n <- nrow(X)
  nK <- length(Ks)
  zc <- stats::qnorm(1 - (1 - conf.level) / 2)

  grid <- expand.grid(profile = seq_len(n), K = Ks, effect = .effect_names_ind,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  effects <- data.frame(profile = grid$profile, K = grid$K,
                        effect = factor(grid$effect, levels = .effect_names_ind),
                        estimate = as.vector(res$est), std.error = as.vector(res$se))
  effects$conf.low <- effects$estimate - zc * effects$std.error
  effects$conf.high <- effects$estimate + zc * effects$std.error
  effects <- effects[order(effects$profile, effects$K, effects$effect), , drop = FALSE]
  rownames(effects) <- NULL

  d <- res$est[, , 6L]
  s <- res$se[, , 6L]
  z <- d / s
  test <- data.frame(profile = rep(seq_len(n), times = nK), K = rep(Ks, each = n),
                     estimate = as.vector(d), std.error = as.vector(s),
                     statistic = as.vector(z), p.value = .p_value(as.vector(z), alternative))
  test <- test[order(test$profile, test$K), , drop = FALSE]
  rownames(test) <- NULL

  gr <- expand.grid(profile = seq_len(n), K = Ks, arm = 1:4, KEEP.OUT.ATTRS = FALSE)
  gor_df <- data.frame(profile = gr$profile, K = gr$K, a1 = .arms[gr$arm, 1L],
                       a2 = .arms[gr$arm, 2L], gor = exp(as.vector(res$lgor)))
  gor_df <- gor_df[order(gor_df$profile, gor_df$K), , drop = FALSE]
  rownames(gor_df) <- NULL

  structure(list(
    effects = effects, test = test, gor = gor_df,
    newdata = if (is.null(newdata)) NULL else as.data.frame(newdata),
    profiles = .profile_labels(object, if (is.null(newdata)) NULL else as.data.frame(newdata), n),
    K = Ks, y1 = object$levels[[1L]][c1 + 1L], y2 = object$levels[[2L]][c2 + 1L],
    response = object$response, alternative = alternative, conf.level = conf.level
  ), class = "biscal_aste")
}

#' Population-level separable treatment effects and tests
#'
#' Estimates the expected absolute separable treatment effects
#' \eqn{EASTE_{j,K} = E\{|STE_{j,K}(X)|\}}, their difference
#' \eqn{DEASTE_K = EASTE_{1,K} - EASTE_{2,K}} and the average effects
#' \eqn{E\{STE_{j,K}(X)\}} and \eqn{E\{TTE_K(X)\}}, where the expectation is
#' over a distribution of the covariates, with delta-method standard errors
#' and a Wald test (Section 3.2.2 of the paper).
#'
#' @section Covariate distribution:
#' The expectation is approximated by the average over the rows of `newdata`.
#' \itemize{
#'   \item By default (`newdata = NULL`) these are the covariates of the
#'     subjects in the data (their empirical distribution). The standard errors
#'     then add the sampling variability of the covariates,
#'     \eqn{\widehat{var}\{h(X)\}/n}, to the delta-method variance of the
#'     estimated coefficients.
#'   \item Alternatively, `newdata` may hold a large number of draws from a
#'     model for the covariate distribution (in the paper, a Gaussian mixture
#'     for the continuous covariates within groups defined by the discrete
#'     ones), or the profiles of a target population. The standard errors then
#'     reflect the uncertainty in the estimated coefficients only, treating
#'     the covariate distribution as known, as in the authors' implementation;
#'     the column `mc_se` gives the Monte Carlo standard error of the average
#'     when the rows are independent draws.
#' }
#'
#' @section Inference:
#' Standard errors use the delta method (the variance formula following
#' equation (27) of the paper), with central-difference derivatives of the
#' average with respect to the coefficients, as in [aste()]. When the average
#' is over the covariates of the subjects in the data (`newdata = NULL`), the
#' sampling variance of the covariates is added (see `covariate_variance`).
#' The paper models the covariate distribution parametrically,
#' \eqn{f_X(x; \alpha)}, and its variance formula also includes the uncertainty
#' of the estimated \eqn{\alpha}; when `newdata` holds draws from such a fitted
#' model, that component is not included, as in the authors' implementation.
#' The default test is that of the paper, \eqn{H_0: DEASTE_K \ge 0} against
#' \eqn{H_1: DEASTE_K < 0}, with p-value \eqn{\Phi(z)}; confidence intervals
#' are Wald intervals.
#'
#' @inheritParams aste
#' @param newdata Optional data frame of covariate values over which the
#'   effects are averaged (see the section on the covariate distribution).
#' @param covariate_variance Add the sampling variance of the covariates to the
#'   standard errors? Defaults to `TRUE` when the empirical distribution of
#'   the data is used (`newdata = NULL`) and `FALSE` otherwise.
#'
#' @return An object of class `"biscal_easte"`, a list with components
#'   \describe{
#'     \item{`effects`}{Data frame with one row per end time and effect
#'       (`ESTE1`, `ESTE2`, `ETTE`, `EASTE1`, `EASTE2`, `DEASTE`): `estimate`,
#'       `std.error`, `conf.low`, `conf.high` and `mc_se`.}
#'     \item{`test`}{Data frame with the test of \eqn{DEASTE_K}: `estimate`,
#'       `std.error`, `statistic` (\eqn{z}) and `p.value`.}
#'     \item{`n_rows`}{Number of covariate rows averaged.}
#'   }
#'   and the settings used.
#'
#' @template ref
#' @seealso [aste()] for effects given covariate values.
#'
#' @examples
#' fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
#'               time = "visit", treatment = "trt")
#' easte(fit)
#'
#' # averaging over draws from a model for the covariates instead
#' set.seed(1)
#' draws <- data.frame(x1 = rnorm(2000), x2 = rbinom(2000, 1, 0.4))
#' easte(fit, newdata = draws)
#'
#' @export
easte <- function(object, newdata = NULL, y1 = NULL, y2 = NULL, K = NULL,
                  alternative = c("less", "greater", "two.sided"), conf.level = 0.95,
                  covariate_variance = is.null(newdata)) {
  .check_biscal(object)
  alternative <- match.arg(alternative)
  .check_conf_level(conf.level)
  .check_flag(covariate_variance, "covariate_variance")
  X <- if (is.null(newdata)) object$X else .newdata_X(object, newdata)
  Ks <- .resolve_K(K, object)
  c1 <- .resolve_threshold(y1, object$levels[[1L]], "y1", object$response[1L])
  c2 <- .resolve_threshold(y2, object$levels[[2L]], "y2", object$response[2L])
  res <- .effects_population(object, X, Ks, c1, c2)
  M <- nrow(X)
  nK <- length(Ks)
  V <- object$vcov
  zc <- stats::qnorm(1 - (1 - conf.level) / 2)

  var_theta <- matrix(NA_real_, nK, 6L)
  for (e in 1:6) {
    G <- matrix(res$grad[, e, ], nrow = nK)
    var_theta[, e] <- pmax(rowSums((G %*% V) * G), 0)
  }
  var_rows <- res$var_rows / M
  se <- sqrt(var_theta + if (covariate_variance) var_rows else 0)
  effects <- data.frame(K = rep(Ks, times = 6L),
                        effect = factor(rep(.effect_names_pop, each = nK),
                                        levels = .effect_names_pop),
                        estimate = as.vector(res$est), std.error = as.vector(se))
  effects$conf.low <- effects$estimate - zc * effects$std.error
  effects$conf.high <- effects$estimate + zc * effects$std.error
  effects$mc_se <- as.vector(sqrt(var_rows))
  effects <- effects[order(effects$K, effects$effect), , drop = FALSE]
  rownames(effects) <- NULL

  z <- res$est[, 6L] / se[, 6L]
  test <- data.frame(K = Ks, estimate = res$est[, 6L], std.error = se[, 6L],
                     statistic = z, p.value = .p_value(z, alternative))

  structure(list(
    effects = effects, test = test, n_rows = M,
    newdata_supplied = !is.null(newdata), covariate_variance = covariate_variance,
    K = Ks, y1 = object$levels[[1L]][c1 + 1L], y2 = object$levels[[2L]][c2 + 1L],
    response = object$response, alternative = alternative, conf.level = conf.level
  ), class = "biscal_easte")
}
