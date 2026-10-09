# Print and convert separable-effect estimates

Methods for the results of
[`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) and
[`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md).

## Usage

``` r
# S3 method for class 'biscal_aste'
print(x, digits = 4, max_print = 6L, ...)

# S3 method for class 'biscal_easte'
print(x, digits = 4, ...)

# S3 method for class 'biscal_aste'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("effects", "test", "gor"),
  ...
)

# S3 method for class 'biscal_easte'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("effects", "test"),
  ...
)
```

## Arguments

- x:

  A result of
  [`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) or
  [`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md).

- digits:

  Number of significant digits.

- max_print:

  Largest number of (profile, end time) combinations printed in full;
  beyond it a compact table is printed.

- ...:

  Further arguments (unused).

- row.names, optional:

  Unused; for compatibility with the generic.

- what:

  `"effects"` (default) for the table of effects, `"test"` for the test
  of the difference of the absolute effects, or `"gor"` (for
  [`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md)
  results) for the generalized odds ratios.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) the
requested table.

## Examples

``` r
fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
              time = "visit", treatment = "trt")
res <- aste(fit, newdata = data.frame(x1 = 0, x2 = 1))
print(res, digits = 3)
#> Separable treatment effects given covariate values (individual level)
#> Thresholds: y1 <= 0, y2 <= 0; Wald 95% confidence intervals
#> 
#> Profile 1 (x1 = 0, x2 = 1), end time K = 10:
#>       Estimate Std. Error Lower 95% Upper 95%
#> STE1   -0.0493     0.0255   -0.0992  0.000648
#> STE2   -0.0485     0.0198   -0.0873  -0.00977
#> TTE    -0.0489     0.0122   -0.0728   -0.0250
#> ASTE1   0.0493     0.0255 -0.000648    0.0992
#> ASTE2   0.0485     0.0198   0.00977    0.0873
#> DASTE 0.000743     0.0385   -0.0747    0.0762
#> H0: DASTE >= 0 vs H1: DASTE < 0: z = 0.0193, p-value = 0.508
as.data.frame(res, what = "test")
#>   profile  K     estimate  std.error  statistic   p.value
#> 1       1 10 0.0007429878 0.03851081 0.01929297 0.5076963
```
