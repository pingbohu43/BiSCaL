# Internal utilities: error helpers, argument checks, a local RNG seed and
# small numerical helpers (softmax, dummy coding, categorical sampling).

.stopf <- function(...) stop(sprintf(...), call. = FALSE)
.warnf <- function(...) warning(sprintf(...), call. = FALSE)

.is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x)
}

.check_number <- function(x, name, positive = FALSE) {
  if (!.is_number(x)) .stopf("`%s` must be a single finite number.", name)
  if (positive && x <= 0) .stopf("`%s` must be positive.", name)
  invisible(as.numeric(x))
}

.check_count <- function(x, name, min = 1L) {
  if (!.is_number(x) || abs(x - round(x)) > 1e-8 || x < min) {
    .stopf("`%s` must be a whole number >= %d.", name, as.integer(min))
  }
  invisible(as.integer(round(x)))
}

.check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    .stopf("`%s` must be TRUE or FALSE.", name)
  }
  invisible(x)
}

.check_column <- function(x, name, data) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    .stopf("`%s` must be the name of a column of `data`, given as a string.", name)
  }
  if (!x %in% names(data)) .stopf("Column '%s' (argument `%s`) is not in `data`.", x, name)
  invisible(x)
}

.check_conf_level <- function(x) {
  if (!.is_number(x) || x <= 0 || x >= 1) .stopf("`conf.level` must be a number in (0, 1).")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Local RNG seed. When `seed` is not NULL the seed is set and the caller's
# random number stream is restored when the calling function exits, so that
# functions of this package never change the user's global RNG state.
# ---------------------------------------------------------------------------
.get_rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

.restore_rng_state <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}

.check_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (!.is_number(seed) || abs(seed - round(seed)) > 0) {
    .stopf("`seed` must be NULL or a single whole number.")
  }
  as.integer(seed)
}

# ---------------------------------------------------------------------------
# Baseline-category logit probabilities.
# `eta` is an n x m matrix of linear predictors for categories 1, ..., m
# (category 0 is the reference with linear predictor 0). Returns the n x (m + 1)
# matrix of probabilities of categories 0, ..., m, computed stably.
# ---------------------------------------------------------------------------
.softmax0 <- function(eta) {
  eta <- as.matrix(eta)
  if (ncol(eta) == 1L) {
    return(cbind(stats::plogis(-eta[, 1L]), stats::plogis(eta[, 1L]), deparse.level = 0))
  }
  mx <- pmax(0, apply(eta, 1L, max))
  e0 <- exp(-mx)
  E <- exp(eta - mx)
  den <- e0 + rowSums(E)
  cbind(e0, E, deparse.level = 0) / den
}

# Dummy coding D(y): n x (L - 1) matrix whose column d equals 1 when y == d
# (y coded 0, ..., L - 1), as in the regression models (12)-(14) of the paper.
.dummy <- function(codes, L) {
  out <- matrix(0, nrow = length(codes), ncol = L - 1L)
  hit <- which(codes > 0L)
  if (length(hit)) out[cbind(hit, codes[hit])] <- 1
  out
}

# Draw one categorical value (0, ..., L - 1) per row of the n x L probability
# matrix P, using one uniform random number per row.
.rcat <- function(P) {
  n <- nrow(P)
  L <- ncol(P)
  u <- stats::runif(n)
  if (L == 2L) return(as.integer(u > P[, 1L]))
  cum <- P[, 1L]
  y <- integer(n)
  for (l in seq_len(L - 1L)) {
    y <- y + as.integer(u > cum)
    if (l < L - 1L) cum <- cum + P[, l + 1L]
  }
  y
}

.format_num <- function(x, digits = 4) {
  out <- formatC(signif(x, digits), digits = digits, format = "g", flag = "#")
  out <- sub("\\.$", "", out)
  out[is.na(x)] <- "NA"
  out
}

.format_p <- function(p) {
  out <- format.pval(p, digits = 3, eps = 1e-16)
  out[is.na(p)] <- "NA"
  out
}
