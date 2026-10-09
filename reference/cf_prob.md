# Counterfactual joint distribution of the two responses

Computes \\P\\Y\_{1,k}(a_1, a_2) = y_1, Y\_{2,k}(a_1, a_2) = y_2 \mid X
= x\\\\, the joint distribution of the two responses at visit \\k\\ had
the treatment components been set to \\(a_1, a_2)\\, through the
recursive formula of Theorems 2 and 3 of the paper, from the fitted (or
given) models. The combinations \\(1, 0)\\ and \\(0, 1)\\ are
counterfactual: they are identified under the conditions of the paper.

## Usage

``` r
cf_prob(object, newdata = NULL, a1, a2, K = NULL)
```

## Arguments

- object:

  A fitted model returned by
  [`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md).

- newdata:

  A data frame of covariate values, one row per covariate profile \\x\\.
  May be omitted when the model has no covariates.

- a1, a2:

  Values (0 or 1) of the treatment components acting on the first and
  the second response.

- K:

  Visit(s) \\k\\ at which the distribution is computed; defaults to the
  last visit of the data. `0` gives the baseline distribution. Visits
  after the last visit of the data are extrapolated with the fitted
  transition models (with a warning).

## Value

An array with dimensions (profile, category of the first response,
category of the second response, visit), with dimnames.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

[`gor()`](https://pingbohu43.github.io/BiSCaL/reference/gor.md) for the
generalized odds ratios.

## Examples

``` r
fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
              time = "visit", treatment = "trt")
pr <- cf_prob(fit, newdata = data.frame(x1 = 0, x2 = 1), a1 = 1, a2 = 0)
pr[1, , , 1]          # joint distribution at the last visit
#>    y2
#> y1           0          1
#>   0 0.02200795 0.03601437
#>   1 0.14134513 0.80063255
sum(pr[1, , , 1])     # sums to one
#> [1] 1
```
