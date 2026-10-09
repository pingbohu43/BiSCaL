# Separable treatment effects and tests given covariate values

Estimates, for given covariate values \\x\\, the separable treatment
effects \\STE\_{1,K}\\ and \\STE\_{2,K}\\ on the two responses, the
total treatment effect \\TTE_K\\, the absolute separable treatment
effects \\ASTE\_{j,K} = \|STE\_{j,K}\|\\ and their difference \\DASTE_K
= ASTE\_{1,K} - ASTE\_{2,K}\\, with delta-method standard errors, and
tests whether the treatment affects the second response more than the
first (Section 3.2.1 of the paper).

## Usage

``` r
aste(
  object,
  newdata = NULL,
  y1 = NULL,
  y2 = NULL,
  K = NULL,
  alternative = c("less", "greater", "two.sided"),
  conf.level = 0.95
)
```

## Arguments

- object:

  A fitted model returned by
  [`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md).

- newdata:

  A data frame of covariate values, one row per covariate profile \\x\\
  (for example, a row of the data, or the smallest or largest value of
  each covariate as in Table 2 of the paper). May be omitted when the
  model has no covariates.

- y1, y2:

  Thresholds of the generalized odds ratio: categories of the first and
  second response (any category except the highest). Default to the
  lowest category, the only choice for binary responses.

- K:

  End time point(s): the visit(s) at which the effects are defined.
  Defaults to the last visit of the data; a vector gives effects at
  several visits (see
  [`plot.biscal_aste()`](https://pingbohu43.github.io/BiSCaL/reference/plot.biscal_aste.md)).
  Visits after the last visit of the data are extrapolated with the
  fitted transition models (with a warning).

- alternative:

  Alternative hypothesis for the test of \\DASTE_K\\ (\\DEASTE_K\\ for
  [`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md)):
  `"less"` (default, as in the paper), `"greater"` or `"two.sided"`.

- conf.level:

  Confidence level of the Wald intervals.

## Value

An object of class `"biscal_aste"`, a list with components

- `effects`:

  Data frame with one row per profile, end time and effect (`STE1`,
  `STE2`, `TTE`, `ASTE1`, `ASTE2`, `DASTE`): `estimate`, `std.error`,
  `conf.low`, `conf.high`.

- `test`:

  Data frame with the test of \\DASTE_K\\ for each profile and end time:
  `estimate`, `std.error`, `statistic` (\\z\\) and `p.value`.

- `gor`:

  Data frame of the estimated generalized odds ratios of the four arms.

and the settings used.

## Definitions

With \\GOR_K(a_1, a_2)\\ the generalized odds ratio of the responses at
the end time point \\K\\ under treatment components \\(a_1, a_2)\\ (see
[`gor()`](https://pingbohu43.github.io/BiSCaL/reference/gor.md)),
\$\$TTE_K = \log\frac{GOR_K(1,1)}{GOR_K(0,0)},\$\$ \$\$STE\_{1,K} =
\log\frac{GOR_K(1,1)}{GOR_K(0,1)} +
\log\frac{GOR_K(1,0)}{GOR_K(0,0)},\$\$ \$\$STE\_{2,K} =
\log\frac{GOR_K(1,1)}{GOR_K(1,0)} +
\log\frac{GOR_K(0,1)}{GOR_K(0,0)},\$\$ so that \\STE\_{1,K} + STE\_{2,K}
= 2\\TTE_K\\ (Theorem 1 of the paper). \\STE\_{1,K}\\ changes the
treatment component acting on the first response with the other
component held fixed, and \\STE\_{2,K}\\ conversely.

## Inference

Standard errors use the delta method with the covariance matrix of the
maximum likelihood estimators (equation (24) of the paper); the
derivatives are central differences with step \\\epsilon^{1/3}\\
(\\\epsilon\\ the machine precision), as in the authors' implementation,
and the absolute values are differentiated through
\\\mathrm{sign}(STE\_{j,K})\\, a subgradient. The default test is that
of the paper, \\H_0: DASTE_K \ge 0\\ against \\H_1: DASTE_K \< 0\\ (the
absolute separable effect on the second response is larger), with
p-value \\\Phi(z)\\, \\z = \widehat{DASTE}\_K / \widehat{se}\\.
Confidence intervals are Wald intervals.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

[`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md) for
the population-level effects;
[`gor()`](https://pingbohu43.github.io/BiSCaL/reference/gor.md);
[`plot.biscal_aste()`](https://pingbohu43.github.io/BiSCaL/reference/plot.biscal_aste.md).

## Examples

``` r
fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
              time = "visit", treatment = "trt")
# effects for two covariate profiles
profiles <- data.frame(x1 = c(-1, 1), x2 = c(0, 1))
res <- aste(fit, newdata = profiles)
res
#> Separable treatment effects given covariate values (individual level)
#> Thresholds: y1 <= 0, y2 <= 0; Wald 95% confidence intervals
#> 
#> Profile 1 (x1 = -1, x2 = 0), end time K = 10:
#>       Estimate Std. Error Lower 95% Upper 95%
#> STE1  -0.07897    0.02238   -0.1228  -0.03509
#> STE2  -0.03184    0.01304  -0.05739 -0.006276
#> TTE   -0.05540    0.01037  -0.07572  -0.03508
#> ASTE1  0.07897    0.02238   0.03509    0.1228
#> ASTE2  0.03184    0.01304  0.006276   0.05739
#> DASTE  0.04713    0.03021  -0.01207    0.1063
#> H0: DASTE >= 0 vs H1: DASTE < 0: z = 1.560, p-value = 0.941
#> 
#> Profile 2 (x1 = 1, x2 = 1), end time K = 10:
#>       Estimate Std. Error Lower 95% Upper 95%
#> STE1  -0.06218    0.02425   -0.1097  -0.01466
#> STE2  -0.05422    0.02212  -0.09758  -0.01086
#> TTE   -0.05820    0.01221  -0.08214  -0.03426
#> ASTE1  0.06218    0.02425   0.01466    0.1097
#> ASTE2  0.05422    0.02212   0.01086   0.09758
#> DASTE 0.007958    0.03947  -0.06940   0.08532
#> H0: DASTE >= 0 vs H1: DASTE < 0: z = 0.2016, p-value = 0.58
as.data.frame(res)
#>    profile  K effect     estimate  std.error     conf.low    conf.high
#> 1        1 10   STE1 -0.078965086 0.02238491 -0.122838708 -0.035091465
#> 2        1 10   STE2 -0.031835515 0.01304067 -0.057394767 -0.006276263
#> 3        1 10    TTE -0.055400301 0.01036754 -0.075720298 -0.035080304
#> 4        1 10  ASTE1  0.078965086 0.02238491  0.035091465  0.122838708
#> 5        1 10  ASTE2  0.031835515 0.01304067  0.006276263  0.057394767
#> 6        1 10  DASTE  0.047129571 0.03020503 -0.012071197  0.106330339
#> 7        2 10   STE1 -0.062179449 0.02424561 -0.109699979 -0.014658919
#> 8        2 10   STE2 -0.054221213 0.02212433 -0.097584108 -0.010858317
#> 9        2 10    TTE -0.058200331 0.01221442 -0.082140161 -0.034260501
#> 10       2 10  ASTE1  0.062179449 0.02424561  0.014658919  0.109699979
#> 11       2 10  ASTE2  0.054221213 0.02212433  0.010858317  0.097584108
#> 12       2 10  DASTE  0.007958236 0.03947028 -0.069402098  0.085318571
```
