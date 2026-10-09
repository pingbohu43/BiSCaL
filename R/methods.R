#' Methods for fitted biscal models
#'
#' Print, summarize and extract information from a model fitted by
#' [biscal()].
#'
#' @param object,x A fitted model returned by [biscal()] (or its summary).
#' @param digits Number of significant digits to print.
#' @param signif.stars Show significance stars in the coefficient tables?
#' @param parm,level Coefficients (names or indices) and confidence level for
#'   [confint()].
#' @param ... Further arguments (unused).
#'
#' @return `coef()` returns the stacked coefficient vector
#'   \eqn{\theta = (\beta_1^\top, \beta_2^\top, \gamma_1^\top, \gamma_2^\top)^\top},
#'   with names `model:term` (for responses with more than two categories
#'   `model[category]:term`); `vcov()` its estimated covariance matrix;
#'   `confint()` Wald confidence intervals; `logLik()` the maximized
#'   log-likelihood; `nobs()` the number of subjects. `summary()` returns an
#'   object of class `"summary.biscal"` whose `coefficients` element is a list of
#'   coefficient tables (estimate, standard error, z value and two-sided
#'   p-value, as in Table 1 of the paper), one per model.
#'
#' @examples
#' fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
#'               time = "visit", treatment = "trt")
#' summary(fit)
#' head(coef(fit))
#' confint(fit, parm = "beta1:trt")
#'
#' @name biscal-methods
NULL

#' @rdname biscal-methods
#' @export
coef.biscal <- function(object, ...) {
  stats::setNames(.flatten(object$coefficients), rownames(object$vcov))
}

#' @rdname biscal-methods
#' @export
vcov.biscal <- function(object, ...) object$vcov

#' @rdname biscal-methods
#' @export
nobs.biscal <- function(object, ...) object$n

#' @rdname biscal-methods
#' @export
logLik.biscal <- function(object, ...) {
  structure(object$loglik, df = nrow(object$vcov), nobs = object$n, class = "logLik")
}

#' @rdname biscal-methods
#' @export
confint.biscal <- function(object, parm, level = 0.95, ...) {
  .check_conf_level(level)
  cf <- coef(object)
  if (missing(parm)) parm <- names(cf)
  if (is.numeric(parm)) parm <- names(cf)[parm]
  bad <- setdiff(parm, names(cf))
  if (length(bad)) .stopf("Unknown coefficient(s): %s.", paste(bad, collapse = ", "))
  se <- sqrt(diag(object$vcov))[parm]
  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  out <- cbind(cf[parm] - z * se, cf[parm] + z * se)
  dimnames(out) <- list(parm, paste(format(100 * c(a, 1 - a), trim = TRUE,
                                           scientific = FALSE, digits = 3), "%"))
  out
}

.model_titles <- function(object) {
  r <- object$response
  c(beta1 = sprintf("Transition model for %s (model (12), j = 1)", r[1L]),
    beta2 = sprintf("Transition model for %s (model (12), j = 2)", r[2L]),
    gamma1 = sprintf("Baseline model for %s (model (13))", r[1L]),
    gamma2 = sprintf("Baseline model for %s (model (14))", r[2L]))
}

.data_lines <- function(object) {
  ntr <- sum(object$A)
  c(sprintf("Subjects: %d (treated %d, control %d); visits 0 to %d (K = %d)",
            object$n, ntr, object$n - ntr, object$K, object$K),
    sprintf("Treatment: %s (control = %s, treated = %s)", object$treatment,
            object$treatment_levels[1L], object$treatment_levels[2L]),
    sprintf("Responses: %s [%s] then %s [%s]", object$response[1L],
            paste(object$levels[[1L]], collapse = " < "), object$response[2L],
            paste(object$levels[[2L]], collapse = " < ")),
    sprintf("Follow-up: mean last visit %.2f; missing after baseline: %s %.1f%%, %s %.1f%%",
            mean(object$last_time), object$response[1L], 100 * object$missing_rate[1L],
            object$response[2L], 100 * object$missing_rate[2L]),
    sprintf("Observations used: transitions %s %d, %s %d; baseline %s %d, %s %d",
            object$response[1L], as.integer(object$nobs_models[["beta1"]]),
            object$response[2L], as.integer(object$nobs_models[["beta2"]]),
            object$response[1L], as.integer(object$nobs_models[["gamma1"]]),
            object$response[2L], as.integer(object$nobs_models[["gamma2"]])))
}

#' @rdname biscal-methods
#' @export
print.biscal <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Separable treatment effects for longitudinal bivariate responses (BiSCaL)\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(paste0(.data_lines(x), collapse = "\n"), "\n", sep = "")
  cat(sprintf("Log-likelihood: %s (df = %d)\n", format(x$loglik, digits = digits + 3L),
              nrow(x$vcov)))
  titles <- .model_titles(x)
  cat("\nCoefficients:\n")
  for (nm in names(x$coefficients)) {
    cat(titles[[nm]], ":\n", sep = "")
    B <- x$coefficients[[nm]]
    if (ncol(B) == 1L) {
      print(stats::setNames(B[, 1L], rownames(B)), digits = digits)
    } else {
      colnames(B) <- paste0("category ", colnames(B))
      print(t(B), digits = digits)
    }
  }
  if (!all(x$converged)) {
    cat("\nWarning: the fit of", paste(names(x$converged)[!x$converged], collapse = ", "),
        "did not converge.\n")
  }
  invisible(x)
}

#' @rdname biscal-methods
#' @export
summary.biscal <- function(object, ...) {
  se_all <- sqrt(diag(object$vcov))
  pos <- 0L
  tabs <- list()
  for (nm in names(object$coefficients)) {
    B <- object$coefficients[[nm]]
    for (d in seq_len(ncol(B))) {
      idx <- pos + seq_len(nrow(B))
      est <- B[, d]
      se <- se_all[idx]
      z <- est / se
      tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = z,
                   `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
      rownames(tab) <- rownames(B)
      key <- if (ncol(B) == 1L) nm else paste0(nm, "[", colnames(B)[d], "]")
      tabs[[key]] <- tab
      pos <- pos + nrow(B)
    }
  }
  structure(list(call = object$call, coefficients = tabs, titles = .model_titles(object),
                 data_lines = .data_lines(object), loglik = object$loglik,
                 df = nrow(object$vcov), levels = object$levels,
                 response = object$response, converged = object$converged),
            class = "summary.biscal")
}

#' @rdname biscal-methods
#' @export
print.summary.biscal <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"), ...) {
  cat("Separable treatment effects for longitudinal bivariate responses (BiSCaL)\n\n")
  cat("Call:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat(paste0(x$data_lines, collapse = "\n"), "\n", sep = "")
  cat(sprintf("Log-likelihood: %s (df = %d)\n", format(x$loglik, digits = digits + 3L), x$df))
  for (key in names(x$coefficients)) {
    nm <- sub("\\[.*$", "", key)
    cat("\n", x$titles[[nm]], sep = "")
    if (grepl("\\[", key)) {
      cat(", category ", sub("^.*\\[(.*)\\]$", "\\1", key), " vs ",
          x$levels[[if (nm %in% c("beta1", "gamma1")) 1L else 2L]][1L], sep = "")
    }
    cat(":\n")
    stats::printCoefmat(x$coefficients[[key]], digits = digits,
                        signif.stars = signif.stars, P.values = TRUE, has.Pvalue = TRUE)
  }
  if (!all(x$converged)) {
    cat("\nWarning: the fit of", paste(names(x$converged)[!x$converged], collapse = ", "),
        "did not converge.\n")
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Effects objects
# ---------------------------------------------------------------------------

.hypothesis_text <- function(alternative, what) {
  switch(alternative,
         less = sprintf("H0: %s >= 0 vs H1: %s < 0", what, what),
         greater = sprintf("H0: %s <= 0 vs H1: %s > 0", what, what),
         two.sided = sprintf("H0: %s = 0 vs H1: %s != 0", what, what))
}

.effects_table <- function(df, conf.level, digits) {
  tab <- data.frame(Estimate = .format_num(df$estimate, digits),
                    `Std. Error` = .format_num(df$std.error, digits),
                    Lower = .format_num(df$conf.low, digits),
                    Upper = .format_num(df$conf.high, digits),
                    check.names = FALSE)
  names(tab)[3:4] <- paste0(c("Lower ", "Upper "), format(100 * conf.level), "%")
  rownames(tab) <- as.character(df$effect)
  tab
}

#' Print and convert separable-effect estimates
#'
#' Methods for the results of [aste()] and [easte()].
#'
#' @param x A result of [aste()] or [easte()].
#' @param digits Number of significant digits.
#' @param max_print Largest number of (profile, end time) combinations printed
#'   in full; beyond it a compact table is printed.
#' @param row.names,optional Unused; for compatibility with the generic.
#' @param what `"effects"` (default) for the table of effects, `"test"` for the
#'   test of the difference of the absolute effects, or `"gor"` (for [aste()]
#'   results) for the generalized odds ratios.
#' @param ... Further arguments (unused).
#'
#' @return `print()` returns `x` invisibly; `as.data.frame()` the requested
#'   table.
#'
#' @examples
#' fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
#'               time = "visit", treatment = "trt")
#' res <- aste(fit, newdata = data.frame(x1 = 0, x2 = 1))
#' print(res, digits = 3)
#' as.data.frame(res, what = "test")
#'
#' @name biscal-effects
NULL

#' @rdname biscal-effects
#' @export
print.biscal_aste <- function(x, digits = 4, max_print = 6L, ...) {
  cat("Separable treatment effects given covariate values (individual level)\n")
  cat(sprintf("Thresholds: %s <= %s, %s <= %s; Wald %s%% confidence intervals\n",
              x$response[1L], x$y1, x$response[2L], x$y2, format(100 * x$conf.level)))
  combos <- unique(x$test[, c("profile", "K")])
  if (nrow(combos) <= max_print) {
    for (i in seq_len(nrow(combos))) {
      pr <- combos$profile[i]
      k <- combos$K[i]
      cat(sprintf("\nProfile %d (%s), end time K = %d:\n", pr, x$profiles[pr], k))
      df <- x$effects[x$effects$profile == pr & x$effects$K == k, , drop = FALSE]
      print(.effects_table(df, x$conf.level, digits))
      tt <- x$test[x$test$profile == pr & x$test$K == k, , drop = FALSE]
      cat(sprintf("%s: z = %s, p-value = %s\n", .hypothesis_text(x$alternative, "DASTE"),
                  .format_num(tt$statistic, digits), .format_p(tt$p.value)))
    }
  } else {
    wide <- stats::reshape(x$effects[x$effects$effect %in% c("STE1", "STE2", "TTE", "DASTE"),
                                     c("profile", "K", "effect", "estimate")],
                           idvar = c("profile", "K"), timevar = "effect", direction = "wide")
    names(wide) <- sub("^estimate\\.", "", names(wide))
    wide <- merge(wide, x$test[, c("profile", "K", "std.error", "p.value")],
                  by = c("profile", "K"))
    wide <- wide[order(wide$profile, wide$K), , drop = FALSE]
    names(wide)[names(wide) == "std.error"] <- "SE(DASTE)"
    for (j in 3:ncol(wide)) wide[[j]] <- signif(wide[[j]], digits)
    cat(sprintf(paste0("\nEstimates of STE1, STE2, TTE and DASTE, with the standard error ",
                       "and p-value of DASTE (%s):\n"),
                .hypothesis_text(x$alternative, "DASTE")))
    print(wide, row.names = FALSE)
    cat("Profiles:\n")
    cat(paste0("  ", seq_along(x$profiles), ": ", x$profiles, collapse = "\n"), "\n")
  }
  invisible(x)
}

#' @rdname biscal-effects
#' @export
print.biscal_easte <- function(x, digits = 4, ...) {
  cat("Separable treatment effects averaged over covariates (population level)\n")
  src <- if (x$newdata_supplied) "the supplied covariate values" else
    "the covariates of the fitted data"
  cat(sprintf("Average over %d rows of %s\n", x$n_rows, src))
  cat(sprintf("Thresholds: %s <= %s, %s <= %s; Wald %s%% confidence intervals\n",
              x$response[1L], x$y1, x$response[2L], x$y2, format(100 * x$conf.level)))
  cat(sprintf("Standard errors %s the sampling variability of the covariates\n",
              if (x$covariate_variance) "include" else "exclude"))
  for (k in x$K) {
    cat(sprintf("\nEnd time K = %d:\n", k))
    df <- x$effects[x$effects$K == k, , drop = FALSE]
    print(.effects_table(df, x$conf.level, digits))
    tt <- x$test[x$test$K == k, , drop = FALSE]
    cat(sprintf("%s: z = %s, p-value = %s\n", .hypothesis_text(x$alternative, "DEASTE"),
                .format_num(tt$statistic, digits), .format_p(tt$p.value)))
  }
  invisible(x)
}

#' @rdname biscal-effects
#' @export
as.data.frame.biscal_aste <- function(x, row.names = NULL, optional = FALSE,
                                      what = c("effects", "test", "gor"), ...) {
  what <- match.arg(what)
  x[[what]]
}

#' @rdname biscal-effects
#' @export
as.data.frame.biscal_easte <- function(x, row.names = NULL, optional = FALSE,
                                       what = c("effects", "test"), ...) {
  what <- match.arg(what)
  x[[what]]
}
