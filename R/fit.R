#' Fit the regression models for separable treatment effects
#'
#' Fits, by maximum likelihood, the four baseline-category logistic regression
#' models through which the separable treatment effects of Hu and Yi (2026) are
#' identified: the transition models for the two responses (model (12) of the
#' paper) and the models for their baseline values (models (13) and (14)).
#' Intermittent missingness and right-censoring are handled as in Section 4 of
#' the paper.
#'
#' @section Model:
#' For subject \eqn{i} at visit \eqn{k = 0, 1, \dots, K}, let \eqn{Y_{1,k}} and
#' \eqn{Y_{2,k}} be the two responses, with \eqn{Y_{1,k}} measured before
#' \eqn{Y_{2,k}}; let \eqn{A} be the treatment and \eqn{X} the baseline
#' covariates. With \eqn{D_j(\cdot)} the dummy (indicator) coding of a
#' previous response, the models are, for every category \eqn{d \ge 1}:
#' \itemize{
#'   \item `beta1` (model (12), \eqn{j = 1}):
#'     \eqn{\log\{P(Y_{1,k} = d \mid \cdot)/P(Y_{1,k} = 0 \mid \cdot)\} =
#'     \beta_{1,0,d} + \beta_{1,1,d}^\top D_1(y_{1,k-1}) +
#'     \beta_{1,2,d}^\top D_2(y_{2,k-1}) + \beta_{1,3,d}\, a + \beta_{1,4,d}^\top x};
#'   \item `beta2` (model (12), \eqn{j = 2}): the same with
#'     \eqn{D_1(y_{1,k})}, the first response at the same visit, in place of
#'     \eqn{D_1(y_{1,k-1})};
#'   \item `gamma1` (model (13)): \eqn{Y_{1,0}} given \eqn{x};
#'   \item `gamma2` (model (14)): \eqn{Y_{2,0}} given \eqn{D_1(y_{1,0})} and
#'     \eqn{x}.
#' }
#' The transition models are common to all visits \eqn{k \ge 1}. Binary
#' responses give ordinary logistic regressions. Every category of each
#' response must occur as an outcome of each model (for example, a response
#' cannot be constant at baseline), and the visits must be numbered
#' consecutively so that transitions between visits \eqn{k - 1} and \eqn{k} are
#' observed.
#'
#' @section Missing data and censoring:
#' Responses that are not observed are coded `NA`, or their rows are omitted
#' (for example, all visits after a subject is censored). Following Section
#' 4.2 of the paper, the transition of the first response at visit \eqn{k}
#' contributes to the likelihood when \eqn{Y_{1,k}}, \eqn{Y_{1,k-1}} and
#' \eqn{Y_{2,k-1}} are observed, and that of the second response when
#' \eqn{Y_{2,k}}, \eqn{Y_{2,k-1}} and \eqn{Y_{1,k}} are observed (the
#' indicators \eqn{m'_{i,j,k}} of the paper). This is valid when missingness
#' and censoring satisfy Conditions AC1 and AC2 of the paper (for example, when
#' they depend only on the treatment, the covariates and earlier responses).
#' The baseline models use the subjects whose baseline responses are observed.
#'
#' @param formula A formula `cbind(y1, y2) ~ x1 + x2 + ...`. The left-hand side
#'   names the two responses, the one measured first at each visit first. They
#'   may be numeric (the sorted distinct values are the ordered categories),
#'   logical, or factors (ordered or not; the level order is used). The
#'   right-hand side lists baseline covariates, which may include factors,
#'   interactions and transformations such as `scale()` or `poly()` (evaluated
#'   on new covariate values with the parameters of the fitted data); use
#'   `~ 1` for none. `offset()` terms are not supported.
#' @param data A data frame in long format: one row per subject and visit.
#' @param id,time,treatment Names (strings) of the columns of `data` holding
#'   the subject identifier, the visit index (`0` for the baseline visit,
#'   `1, ..., K` for the scheduled follow-up visits) and the treatment. The
#'   treatment may be coded 0/1, logical, or a two-level factor whose first
#'   level is the control group. The treatment and covariates must be constant
#'   within subjects.
#' @param control A list of convergence settings: `epsilon` (relative change in
#'   the log-likelihood, default `1e-8`) and `maxit` (maximum number of
#'   iterations, default `25`; at least 100 are allowed for responses with more
#'   than two categories).
#'
#' @return An object of class `"biscal"`, a list with components including
#'   \describe{
#'     \item{`coefficients`}{List of the four coefficient matrices `beta1`,
#'       `beta2`, `gamma1` and `gamma2` (terms in rows, categories
#'       \eqn{d = 1, \dots, L - 1} in columns).}
#'     \item{`vcov`}{Block-diagonal covariance matrix of all coefficients
#'       (inverse observed information), in the order of [coef.biscal()].}
#'     \item{`K`}{The last visit index \eqn{K}, the default end time point of
#'       the effects.}
#'     \item{`n`, `nobs_models`, `loglik`}{Number of subjects, number of
#'       observations used by each model and the maximized log-likelihood.}
#'     \item{`X`, `A`, `covariate_data`}{Subject-level covariate matrix
#'       (used by [easte()] by default), treatment and covariate values.}
#'   }
#'   and information on the data used by the methods.
#'
#' @template ref
#' @seealso [biscal-methods] for `summary()`, `coef()` and other methods;
#'   [aste()] and [easte()] for the separable treatment effects and tests.
#'
#' @examples
#' fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
#'               time = "visit", treatment = "trt")
#' fit
#' summary(fit)
#'
#' @export
biscal <- function(formula, data, id, time, treatment, control = list()) {
  cl <- match.call()
  ctrl <- .biscal_control(control)
  d <- .prepare_data(formula, data, id, time, treatment)
  des <- .build_designs(d)
  r <- d$resp_names
  labels <- c(beta1 = sprintf("transition model for %s", r[1L]),
              beta2 = sprintf("transition model for %s", r[2L]),
              gamma1 = sprintf("baseline model for %s", r[1L]),
              gamma2 = sprintf("baseline model for %s", r[2L]))
  fits <- lapply(names(des), function(nm) {
    ds <- des[[nm]]
    .fit_mlogit(ds$Z, ds$y, ds$L, ds$lev, labels[[nm]], ctrl)
  })
  names(fits) <- names(des)
  coefs <- lapply(names(fits), function(nm) {
    B <- fits[[nm]]$coef
    dimnames(B) <- list(colnames(des[[nm]]$Z), des[[nm]]$lev[-1L])
    B
  })
  names(coefs) <- names(fits)
  theta_names <- .theta_names(coefs)
  V <- .block_diag(lapply(fits, `[[`, "vcov"))
  dimnames(V) <- list(theta_names, theta_names)

  structure(list(
    coefficients = coefs,
    vcov = V,
    call = cl,
    formula = formula,
    terms = d$terms,
    xlevels = d$xlevels,
    contrasts = d$contrasts,
    cov_vars = d$cov_vars,
    response = r,
    levels = list(d$lev1, d$lev2),
    treatment = d$treatment,
    treatment_levels = d$treatment_levels,
    id = d$id,
    time = d$time,
    K = d$K,
    n = d$n,
    ids = d$ids,
    X = d$X,
    A = d$A,
    covariate_data = d$covariate_data,
    last_time = d$last_time,
    observed = list(y1 = colSums(!is.na(d$Y1w)), y2 = colSums(!is.na(d$Y2w))),
    missing_rate = .missing_rates(d),
    nobs_models = vapply(fits, `[[`, numeric(1), "nobs"),
    loglik_models = vapply(fits, `[[`, numeric(1), "loglik"),
    loglik = sum(vapply(fits, `[[`, numeric(1), "loglik")),
    converged = vapply(fits, `[[`, logical(1), "converged"),
    model_labels = labels,
    control = ctrl
  ), class = "biscal")
}

# Proportion of missing values of each response over the visits 1, ..., C_i
# before each subject's last recorded visit C_i.
.missing_rates <- function(d) {
  ok <- d$last_time >= 1L
  if (!any(ok)) return(c(NA_real_, NA_real_))
  tot <- sum(d$last_time[ok])
  within <- outer(seq_len(d$n), seq_len(d$K + 1L), function(i, j) j - 1L <= d$last_time[i]) &
    col(d$Y1w) > 1L
  c(sum(is.na(d$Y1w) & within) / tot, sum(is.na(d$Y2w) & within) / tot)
}

# Names of the stacked coefficient vector theta = (beta1, beta2, gamma1, gamma2),
# each block ordered by outcome category, then by term.
.theta_names <- function(coefs) {
  unlist(lapply(names(coefs), function(nm) {
    B <- coefs[[nm]]
    if (ncol(B) == 1L) {
      paste0(nm, ":", rownames(B))
    } else {
      as.vector(outer(rownames(B), colnames(B), function(term, cat) {
        paste0(nm, "[", cat, "]:", term)
      }))
    }
  }), use.names = FALSE)
}

.block_diag <- function(mats) {
  sizes <- vapply(mats, nrow, integer(1))
  out <- matrix(0, sum(sizes), sum(sizes))
  pos <- 0L
  for (M in mats) {
    i <- pos + seq_len(nrow(M))
    out[i, i] <- M
    pos <- pos + nrow(M)
  }
  out
}

.flatten <- function(coefs) unlist(lapply(coefs, as.vector), use.names = FALSE)

.unflatten <- function(theta, template) {
  out <- template
  pos <- 0L
  for (nm in names(template)) {
    len <- length(template[[nm]])
    out[[nm]][] <- theta[pos + seq_len(len)]
    pos <- pos + len
  }
  out
}
