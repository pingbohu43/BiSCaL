# Plot separable treatment effects

Plots the estimated effects with Wald confidence intervals. When the
effects were computed at several end time points (a vector `K` in
[`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) or
[`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md)),
they are drawn as curves over the end time with pointwise confidence
bands; otherwise as points with interval bars.

## Usage

``` r
# S3 method for class 'biscal_aste'
plot(
  x,
  effects = c("STE1", "STE2", "TTE"),
  profiles = NULL,
  band = TRUE,
  col = NULL,
  main = NULL,
  xlab = "End time K",
  ylab = "Effect (log GOR ratio)",
  ylim = NULL,
  legend = TRUE,
  mfrow = NULL,
  ...
)

# S3 method for class 'biscal_easte'
plot(
  x,
  effects = c("ESTE1", "ESTE2", "ETTE"),
  band = TRUE,
  col = NULL,
  main = "",
  xlab = "End time K",
  ylab = "Average effect (log GOR ratio)",
  ylim = NULL,
  legend = TRUE,
  ...
)
```

## Arguments

- x:

  A result of
  [`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) or
  [`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md).

- effects:

  Effects to draw. For
  [`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md)
  results any of `"STE1"`, `"STE2"`, `"TTE"`, `"ASTE1"`, `"ASTE2"` and
  `"DASTE"`; for
  [`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md)
  results any of `"ESTE1"`, `"ESTE2"`, `"ETTE"`, `"EASTE1"`, `"EASTE2"`
  and `"DEASTE"`.

- profiles:

  For [`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md)
  results, the covariate profiles (row numbers of `newdata`) to draw,
  one panel each.

- band:

  Draw confidence intervals?

- col:

  Colors, one per effect.

- main:

  Panel titles; by default the covariate values of each profile.

- xlab, ylab:

  Axis labels.

- ylim:

  Vertical range; by default it covers the intervals and zero.

- legend:

  Draw a legend (for curves over the end time)?

- mfrow:

  Panel layout `c(rows, columns)` for several profiles.

- ...:

  Further graphical arguments passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`NULL`, invisibly. Called for its side effect.

## Examples

``` r
fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
              time = "visit", treatment = "trt")
res <- aste(fit, newdata = data.frame(x1 = 0, x2 = 1), K = 1:10)
plot(res)

plot(res, effects = "DASTE")

```
