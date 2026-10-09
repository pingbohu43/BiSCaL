# Data preparation: validation of the long-format data, coding of the
# responses and the treatment, subject-level covariates, the subject-by-visit
# response matrices and the four regression data sets of the likelihood.

# Code a response as integers 0, ..., L - 1 and return its category labels.
.code_response <- function(y, name) {
  if (is.factor(y)) {
    used <- levels(droplevels(y))
    if (length(used) < length(levels(y))) {
      .warnf("Unused categories of response '%s' were dropped: %s.", name,
             paste(setdiff(levels(y), used), collapse = ", "))
    }
    lev <- used
    codes <- match(as.character(y), lev) - 1L
  } else if (is.logical(y)) {
    lev <- c("FALSE", "TRUE")
    codes <- as.integer(y)
    if (!any(y %in% TRUE) || !any(y %in% FALSE)) lev <- character(0)
  } else if (is.numeric(y)) {
    u <- sort(unique(y[!is.na(y)]))
    lev <- as.character(u)
    codes <- match(y, u) - 1L
  } else {
    .stopf(paste0("Response '%s' must be numeric, logical or a factor. Convert character ",
                  "responses to an (ordered) factor so that the order of the categories ",
                  "is explicit."), name)
  }
  if (length(lev) < 2L) {
    .stopf("Response '%s' must have at least two observed categories.", name)
  }
  list(codes = as.integer(codes), levels = lev)
}

# Code the treatment as 0 (control) / 1 (treated).
.code_treatment <- function(a, name) {
  if (is.logical(a)) {
    codes <- as.integer(a)
    lev <- c("FALSE", "TRUE")
  } else if (is.factor(a)) {
    lev <- levels(droplevels(a))
    if (length(lev) != 2L) {
      .stopf("Treatment '%s' must have exactly two groups (found %d).", name, length(lev))
    }
    codes <- match(as.character(a), lev) - 1L
  } else if (is.numeric(a)) {
    if (!all(a %in% c(0, 1))) {
      .stopf(paste0("Treatment '%s' must be coded 0 (control) and 1 (treated), or be a ",
                    "logical or a two-level factor."), name)
    }
    codes <- as.integer(a)
    lev <- c("0", "1")
  } else {
    .stopf(paste0("Treatment '%s' must be numeric 0/1, logical or a factor whose first ",
                  "level is the control group."), name)
  }
  list(codes = as.integer(codes), levels = lev)
}

# Is `v` constant within subjects (NA-aware)? `first` = row of the first
# record of each subject, `idx` = subject index of every row.
.constant_within <- function(v, idx, first) {
  if (is.factor(v)) v <- as.character(v)
  fv <- v[first][idx]
  same <- (is.na(v) & is.na(fv)) | (!is.na(v) & !is.na(fv) & v == fv)
  all(same)
}

.prepare_data <- function(formula, data, id, time, treatment) {
  if (!is.data.frame(data)) .stopf("`data` must be a data frame.")
  if (!nrow(data)) .stopf("`data` has no rows.")
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    .stopf("`formula` must be a two-sided formula of the form cbind(y1, y2) ~ covariates.")
  }
  .check_column(id, "id", data)
  .check_column(time, "time", data)
  .check_column(treatment, "treatment", data)
  if (anyDuplicated(c(id, time, treatment))) {
    .stopf("`id`, `time` and `treatment` must name three different columns.")
  }

  # responses -----------------------------------------------------------------
  lhs <- formula[[2L]]
  if (!is.call(lhs) || !identical(lhs[[1L]], as.name("cbind")) || length(lhs) != 3L) {
    .stopf(paste0("The left-hand side of `formula` must be cbind(<first response>, ",
                  "<second response>), with the response measured first at each visit ",
                  "listed first."))
  }
  env <- environment(formula)
  if (is.null(env)) env <- parent.frame()
  resp_names <- c(paste(deparse(lhs[[2L]]), collapse = ""),
                  paste(deparse(lhs[[3L]]), collapse = ""))
  if (resp_names[1L] == resp_names[2L]) .stopf("The two responses must be different.")
  y1 <- eval(lhs[[2L]], data, env)
  y2 <- eval(lhs[[3L]], data, env)
  if (length(y1) != nrow(data) || length(y2) != nrow(data)) {
    .stopf("Each response must have one value per row of `data`.")
  }
  r1 <- .code_response(y1, resp_names[1L])
  r2 <- .code_response(y2, resp_names[2L])

  # covariates (right-hand side) ------------------------------------------------
  rhs <- formula[-2L]
  cov_vars <- all.vars(rhs)
  if ("." %in% cov_vars) {
    .stopf("`.` is not supported in `formula`; list the baseline covariates explicitly.")
  }
  tt <- stats::terms(rhs)
  if (!is.null(attr(tt, "offset"))) .stopf("`offset()` terms are not supported in `formula`.")
  if (attr(tt, "intercept") == 0L) {
    .stopf("All models include an intercept; remove `- 1` or `+ 0` from `formula`.")
  }
  bad <- intersect(cov_vars, c(id, time, treatment))
  if (length(bad)) {
    .stopf(paste0("'%s' cannot be used as a covariate: the identifier, visit time and ",
                  "treatment are specified by `id`, `time` and `treatment`."), bad[1L])
  }
  bad <- intersect(cov_vars, all.vars(lhs))
  if (length(bad)) .stopf("Response '%s' cannot also be a covariate.", bad[1L])
  missing_cols <- setdiff(cov_vars, names(data))
  if (length(missing_cols)) {
    .stopf("Covariate '%s' is not a column of `data`.", missing_cols[1L])
  }

  # id, time, treatment -----------------------------------------------------------
  idv <- data[[id]]
  tv <- data[[time]]
  if (anyNA(idv)) .stopf("Column '%s' (`id`) contains missing values.", id)
  if (anyNA(tv)) .stopf("Column '%s' (`time`) contains missing values.", time)
  if (!is.numeric(tv) || any(abs(tv - round(tv)) > 1e-8) || any(tv < 0)) {
    .stopf(paste0("Column '%s' (`time`) must contain the visit index: 0 for the baseline ",
                  "visit and 1, 2, ..., K for the follow-up visits."), time)
  }
  tv <- as.integer(round(tv))
  if (min(tv) != 0L) {
    .stopf(paste0("No baseline visit found: column '%s' (`time`) must code the baseline ",
                  "visit as 0 and the follow-up visits as 1, 2, ..., K."), time)
  }
  K <- max(tv)
  if (K < 1L) .stopf("At least one follow-up visit (time >= 1) is needed.")
  subj <- factor(as.character(idv), levels = unique(as.character(idv)))
  idx <- as.integer(subj)
  n <- nlevels(subj)
  if (anyDuplicated(cbind(idx, tv))) {
    .stopf("Some subjects have more than one row for the same visit time.")
  }
  first <- match(seq_len(n), idx)

  av <- data[[treatment]]
  if (anyNA(av)) .stopf("Column '%s' (`treatment`) contains missing values.", treatment)
  if (!.constant_within(av, idx, first)) {
    .stopf("The treatment '%s' must be constant within subjects.", treatment)
  }
  tr <- .code_treatment(av, treatment)
  A <- tr$codes[first]
  if (length(unique(A)) < 2L) .stopf("Both treatment groups must be present in the data.")

  for (v in cov_vars) {
    if (!.constant_within(data[[v]], idx, first)) {
      .stopf(paste0("Covariate '%s' varies within subjects; only baseline (time-invariant) ",
                    "covariates are supported."), v)
    }
  }
  sdat <- data[first, cov_vars, drop = FALSE]
  rownames(sdat) <- NULL
  p_terms <- length(attr(tt, "term.labels"))
  if (p_terms == 0L) {
    X <- matrix(numeric(0), nrow = n, ncol = 0L)
    xlevels <- NULL
    contrasts <- NULL
  } else {
    mf <- stats::model.frame(tt, data = sdat, na.action = stats::na.pass,
                             drop.unused.levels = TRUE)
    nas <- !stats::complete.cases(mf)
    if (any(nas)) {
      .stopf(paste0("Covariates are missing for %d subject(s); remove these subjects or ",
                    "impute their covariates before fitting."), sum(nas))
    }
    # the terms of the model frame keep the information ("predvars") needed to
    # evaluate data-dependent transformations such as scale() or poly() on new data
    tt <- attr(mf, "terms")
    X <- stats::model.matrix(tt, mf)
    xlevels <- stats::.getXlevels(tt, mf)
    contrasts <- attr(X, "contrasts")
    if (qr(X)$rank < ncol(X)) {
      .stopf(paste0("The covariates are collinear (or constant), so their coefficients ",
                    "are not identifiable."))
    }
    X <- X[, -1L, drop = FALSE]
    attr(X, "assign") <- NULL
    attr(X, "contrasts") <- NULL
  }
  rownames(X) <- NULL

  # subject-by-visit response matrices (column k + 1 = visit k) ------------------
  Y1w <- matrix(NA_integer_, nrow = n, ncol = K + 1L)
  Y2w <- matrix(NA_integer_, nrow = n, ncol = K + 1L)
  Y1w[cbind(idx, tv + 1L)] <- r1$codes
  Y2w[cbind(idx, tv + 1L)] <- r2$codes
  last_time <- as.integer(tapply(tv, idx, max))

  list(Y1w = Y1w, Y2w = Y2w, A = A, X = X, ids = levels(subj), n = n, K = K,
       lev1 = r1$levels, lev2 = r2$levels, L1 = length(r1$levels), L2 = length(r2$levels),
       resp_names = resp_names, terms = tt, xlevels = xlevels, contrasts = contrasts,
       cov_vars = cov_vars, covariate_data = sdat, treatment_levels = tr$levels,
       last_time = last_time, id = id, time = time, treatment = treatment)
}

# Column of ones (intercept) for the rows `idx`, also when there are none.
.ones <- function(idx) matrix(1, nrow = length(idx), ncol = 1L)

# Regression data sets of the four factors of the likelihood (Section 4.2 of
# the paper). A transition of response 1 at visit k contributes when
# Y1[k], Y1[k-1] and Y2[k-1] are observed (m'_{1,k} = 1); a transition of
# response 2 contributes when Y2[k], Y2[k-1] and Y1[k] are observed
# (m'_{2,k} = 1). Rows are stacked by visit, then by subject.
.build_designs <- function(d) {
  K <- d$K
  L1 <- d$L1
  L2 <- d$L2
  X <- d$X
  A <- d$A
  Y1w <- d$Y1w
  Y2w <- d$Y2w
  r1 <- d$resp_names[1L]
  r2 <- d$resp_names[2L]
  xn <- colnames(X)
  if (is.null(xn)) xn <- character(0)
  nm_b1 <- paste0(r1, "[0]=", d$lev1[-1L])
  nm_lag1 <- paste0(r1, "[k-1]=", d$lev1[-1L])
  nm_cur1 <- paste0(r1, "[k]=", d$lev1[-1L])
  nm_lag2 <- paste0(r2, "[k-1]=", d$lev2[-1L])

  o1 <- which(!is.na(Y1w[, 1L]))
  W1 <- cbind(.ones(o1), X[o1, , drop = FALSE])
  colnames(W1) <- c("(Intercept)", xn)
  o2 <- which(!is.na(Y1w[, 1L]) & !is.na(Y2w[, 1L]))
  W2 <- cbind(.ones(o2), .dummy(Y1w[o2, 1L], L1), X[o2, , drop = FALSE])
  colnames(W2) <- c("(Intercept)", nm_b1, xn)

  rows1 <- lapply(seq_len(K), function(k) {
    which(!is.na(Y1w[, k + 1L]) & !is.na(Y1w[, k]) & !is.na(Y2w[, k]))
  })
  s1 <- unlist(rows1, use.names = FALSE)
  k1 <- rep.int(seq_len(K), lengths(rows1))
  Z1 <- cbind(.ones(s1), .dummy(Y1w[cbind(s1, k1)], L1), .dummy(Y2w[cbind(s1, k1)], L2),
              matrix(A[s1], ncol = 1L), X[s1, , drop = FALSE])
  colnames(Z1) <- c("(Intercept)", nm_lag1, nm_lag2, d$treatment, xn)

  rows2 <- lapply(seq_len(K), function(k) {
    which(!is.na(Y2w[, k + 1L]) & !is.na(Y2w[, k]) & !is.na(Y1w[, k + 1L]))
  })
  s2 <- unlist(rows2, use.names = FALSE)
  k2 <- rep.int(seq_len(K), lengths(rows2))
  Z2 <- cbind(.ones(s2), .dummy(Y1w[cbind(s2, k2 + 1L)], L1), .dummy(Y2w[cbind(s2, k2)], L2),
              matrix(A[s2], ncol = 1L), X[s2, , drop = FALSE])
  colnames(Z2) <- c("(Intercept)", nm_cur1, nm_lag2, d$treatment, xn)

  list(beta1 = list(Z = Z1, y = Y1w[cbind(s1, k1 + 1L)], L = L1, lev = d$lev1),
       beta2 = list(Z = Z2, y = Y2w[cbind(s2, k2 + 1L)], L = L2, lev = d$lev2),
       gamma1 = list(Z = W1, y = Y1w[o1, 1L], L = L1, lev = d$lev1),
       gamma2 = list(Z = W2, y = Y2w[o2, 1L], L = L2, lev = d$lev2))
}
