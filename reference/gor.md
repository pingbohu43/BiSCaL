# Generalized odds ratios of the two responses

Computes the generalized odds ratio \\GOR_k(y_1, y_2, a_1, a_2, x)\\ of
equation (2) of the paper: the odds ratio of the \\2 \times 2\\ table
obtained by dichotomizing the counterfactual responses \\Y\_{1,k}(a_1,
a_2)\\ at \\y_1\\ and \\Y\_{2,k}(a_1, a_2)\\ at \\y_2\\, given \\X =
x\\, \$\$\frac{P(Y_1 \le y_1, Y_2 \le y_2) P(Y_1 \> y_1, Y_2 \>
y_2)}{P(Y_1 \le y_1, Y_2 \> y_2) P(Y_1 \> y_1, Y_2 \le y_2)}.\$\$ Values
above one indicate that the two responses tend to fall on the same side
of their thresholds.

## Usage

``` r
gor(object, newdata = NULL, a1 = 1, a2 = 1, y1 = NULL, y2 = NULL, K = NULL)
```

## Arguments

- object:

  A fitted model returned by
  [`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md).

- newdata:

  A data frame of covariate values, one row per covariate profile \\x\\.
  May be omitted when the model has no covariates.

- a1, a2:

  Values (0 or 1) of the treatment components; vectors of the same
  length give several arms at once.

- y1, y2:

  Thresholds: categories of the first and second response (any category
  except the highest). Default to the lowest category, which is the only
  choice for binary responses.

- K:

  Visit(s) \\k \ge 1\\ at which the ratios are computed; defaults to the
  last visit of the data. Visits after the last visit of the data are
  extrapolated with the fitted transition models (with a warning).

## Value

A data frame with columns `profile`, `K`, `a1`, `a2` and `gor`.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

[`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) for
the separable treatment effects built from these ratios.

## Examples

``` r
fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
              time = "visit", treatment = "trt")
gor(fit, newdata = data.frame(x1 = 0, x2 = 1), a1 = c(1, 1, 0, 0),
    a2 = c(1, 0, 1, 0))
#>   profile  K a1 a2      gor
#> 1       1 10  1  1 3.389511
#> 2       1 10  1  0 3.461436
#> 3       1 10  0  1 3.462722
#> 4       1 10  0  0 3.559347
```
