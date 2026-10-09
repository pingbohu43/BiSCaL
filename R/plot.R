# Graphics (base R).

# Categorical colors in a fixed order (color-vision-deficiency checked).
.biscal_palette <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#4a3aa7")

.col_alpha <- function(col, alpha) {
  rgb <- grDevices::col2rgb(col) / 255
  grDevices::rgb(rgb[1L, ], rgb[2L, ], rgb[3L, ], alpha = alpha)
}

# Draw the effects in `df` (columns K, effect, estimate, conf.low, conf.high):
# lines over the end time with confidence bands when there are several end
# times, otherwise points with interval bars.
.plot_effect_panel <- function(df, effects, col, band, main, xlab, ylab, ylim, legend, ...) {
  df <- df[df$effect %in% effects, , drop = FALSE]
  Ks <- sort(unique(df$K))
  if (is.null(ylim)) {
    rng <- range(c(df$estimate, if (band) c(df$conf.low, df$conf.high), 0), finite = TRUE)
    ylim <- rng
  }
  if (length(Ks) > 1L) {
    graphics::plot(NA, xlim = range(Ks), ylim = ylim, xlab = xlab, ylab = ylab, main = main,
                   ...)
    graphics::abline(h = 0, col = "grey60", lty = 3)
    for (i in seq_along(effects)) {
      d <- df[df$effect == effects[i], , drop = FALSE]
      d <- d[order(d$K), , drop = FALSE]
      if (band) {
        graphics::polygon(c(d$K, rev(d$K)), c(d$conf.low, rev(d$conf.high)),
                          col = .col_alpha(col[i], 0.18), border = NA)
      }
      graphics::lines(d$K, d$estimate, col = col[i], lwd = 2)
    }
  } else {
    pos <- seq_along(effects)
    graphics::plot(NA, xlim = c(0.5, length(effects) + 0.5), ylim = ylim, xaxt = "n",
                   xlab = "", ylab = ylab, main = main, ...)
    graphics::axis(1, at = pos, labels = effects)
    graphics::abline(h = 0, col = "grey60", lty = 3)
    for (i in seq_along(effects)) {
      d <- df[df$effect == effects[i], , drop = FALSE]
      if (band) {
        graphics::segments(pos[i], d$conf.low, pos[i], d$conf.high, col = col[i], lwd = 2)
      }
      graphics::points(pos[i], d$estimate, pch = 19, col = col[i], cex = 1.3)
    }
  }
  if (legend && length(Ks) > 1L) {
    graphics::legend("topright", legend = effects, col = col, lwd = 2, bty = "n",
                     inset = 0.01)
  }
  invisible(NULL)
}

#' Plot separable treatment effects
#'
#' Plots the estimated effects with Wald confidence intervals. When the effects
#' were computed at several end time points (a vector `K` in [aste()] or
#' [easte()]), they are drawn as curves over the end time with pointwise
#' confidence bands; otherwise as points with interval bars.
#'
#' @param x A result of [aste()] or [easte()].
#' @param effects Effects to draw. For [aste()] results any of `"STE1"`,
#'   `"STE2"`, `"TTE"`, `"ASTE1"`, `"ASTE2"` and `"DASTE"`; for [easte()]
#'   results any of `"ESTE1"`, `"ESTE2"`, `"ETTE"`, `"EASTE1"`, `"EASTE2"` and
#'   `"DEASTE"`.
#' @param profiles For [aste()] results, the covariate profiles (row numbers
#'   of `newdata`) to draw, one panel each.
#' @param band Draw confidence intervals?
#' @param col Colors, one per effect.
#' @param main Panel titles; by default the covariate values of each profile.
#' @param xlab,ylab Axis labels.
#' @param ylim Vertical range; by default it covers the intervals and zero.
#' @param legend Draw a legend (for curves over the end time)?
#' @param mfrow Panel layout `c(rows, columns)` for several profiles.
#' @param ... Further graphical arguments passed to [graphics::plot()].
#'
#' @return `NULL`, invisibly. Called for its side effect.
#'
#' @examples
#' fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
#'               time = "visit", treatment = "trt")
#' res <- aste(fit, newdata = data.frame(x1 = 0, x2 = 1), K = 1:10)
#' plot(res)
#' plot(res, effects = "DASTE")
#'
#' @export
plot.biscal_aste <- function(x, effects = c("STE1", "STE2", "TTE"), profiles = NULL,
                             band = TRUE, col = NULL, main = NULL, xlab = "End time K",
                             ylab = "Effect (log GOR ratio)", ylim = NULL, legend = TRUE,
                             mfrow = NULL, ...) {
  effects <- match.arg(effects, .effect_names_ind, several.ok = TRUE)
  .check_flag(band, "band")
  .check_flag(legend, "legend")
  all_pr <- sort(unique(x$effects$profile))
  if (is.null(profiles)) profiles <- all_pr
  if (!all(profiles %in% all_pr)) .stopf("`profiles` must be row numbers of `newdata`.")
  if (is.null(col)) col <- .biscal_palette[match(effects, .effect_names_ind)]
  col <- rep_len(col, length(effects))
  if (is.null(main)) main <- if (length(profiles) == 1L && x$profiles[1L] == "(no covariates)")
    "" else x$profiles[profiles]
  main <- rep_len(main, length(profiles))
  if (length(profiles) > 1L) {
    if (is.null(mfrow)) mfrow <- c(1L, length(profiles))
    op <- graphics::par(mfrow = mfrow)
    on.exit(graphics::par(op), add = TRUE)
  }
  for (i in seq_along(profiles)) {
    df <- x$effects[x$effects$profile == profiles[i], , drop = FALSE]
    .plot_effect_panel(df, effects, col, band, main[i], xlab, ylab, ylim, legend, ...)
  }
  invisible(NULL)
}

#' @rdname plot.biscal_aste
#' @export
plot.biscal_easte <- function(x, effects = c("ESTE1", "ESTE2", "ETTE"), band = TRUE,
                              col = NULL, main = "", xlab = "End time K",
                              ylab = "Average effect (log GOR ratio)", ylim = NULL,
                              legend = TRUE, ...) {
  effects <- match.arg(effects, .effect_names_pop, several.ok = TRUE)
  .check_flag(band, "band")
  .check_flag(legend, "legend")
  if (is.null(col)) col <- .biscal_palette[match(effects, .effect_names_pop)]
  col <- rep_len(col, length(effects))
  .plot_effect_panel(x$effects, effects, col, band, main, xlab, ylab, ylim, legend, ...)
  invisible(NULL)
}
