# Population-level separable treatment effects and tests

Estimates the expected absolute separable treatment effects
\\EASTE\_{j,K} = E\\\|STE\_{j,K}(X)\|\\\\, their difference \\DEASTE_K =
EASTE\_{1,K} - EASTE\_{2,K}\\ and the average effects
\\E\\STE\_{j,K}(X)\\\\ and \\E\\TTE_K(X)\\\\, where the expectation is
over a distribution of the covariates, with delta-method standard errors
and a Wald test (Section 3.2.2 of the paper).

## Usage

``` r
easte(
  object,
  newdata = NULL,
  y1 = NULL,
  y2 = NULL,
  K = NULL,
  alternative = c("less", "greater", "two.sided"),
  conf.level = 0.95,
  covariate_variance = is.null(newdata)
)
```

## Arguments

- object:

  A fitted model returned by
  [`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md).

- newdata:

  Optional data frame of covariate values over which the effects are
  averaged (see the section on the covariate distribution).

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
  `easte()`): `"less"` (default, as in the paper), `"greater"` or
  `"two.sided"`.

- conf.level:

  Confidence level of the Wald intervals.

- covariate_variance:

  Add the sampling variance of the covariates to the standard errors?
  Defaults to `TRUE` when the empirical distribution of the data is used
  (`newdata = NULL`) and `FALSE` otherwise.

## Value

An object of class `"biscal_easte"`, a list with components

- `effects`:

  Data frame with one row per end time and effect (`ESTE1`, `ESTE2`,
  `ETTE`, `EASTE1`, `EASTE2`, `DEASTE`): `estimate`, `std.error`,
  `conf.low`, `conf.high` and `mc_se`.

- `test`:

  Data frame with the test of \\DEASTE_K\\: `estimate`, `std.error`,
  `statistic` (\\z\\) and `p.value`.

- `n_rows`:

  Number of covariate rows averaged.

and the settings used.

## Covariate distribution

The expectation is approximated by the average over the rows of
`newdata`.

- By default (`newdata = NULL`) these are the covariates of the subjects
  in the data (their empirical distribution). The standard errors then
  add the sampling variability of the covariates,
  \\\widehat{var}\\h(X)\\/n\\, to the delta-method variance of the
  estimated coefficients.

- Alternatively, `newdata` may hold a large number of draws from a model
  for the covariate distribution (in the paper, a Gaussian mixture for
  the continuous covariates within groups defined by the discrete ones),
  or the profiles of a target population. The standard errors then
  reflect the uncertainty in the estimated coefficients only, treating
  the covariate distribution as known, as in the authors'
  implementation; the column `mc_se` gives the Monte Carlo standard
  error of the average when the rows are independent draws.

## Inference

Standard errors use the delta method (the variance formula following
equation (27) of the paper), with central-difference derivatives of the
average with respect to the coefficients, as in
[`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md). When
the average is over the covariates of the subjects in the data
(`newdata = NULL`), the sampling variance of the covariates is added
(see `covariate_variance`). The paper models the covariate distribution
parametrically, \\f_X(x; \alpha)\\, and its variance formula also
includes the uncertainty of the estimated \\\alpha\\; when `newdata`
holds draws from such a fitted model, that component is not included, as
in the authors' implementation. The default test is that of the paper,
\\H_0: DEASTE_K \ge 0\\ against \\H_1: DEASTE_K \< 0\\, with p-value
\\\Phi(z)\\; confidence intervals are Wald intervals.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

[`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) for
effects given covariate values.

## Examples

``` r
fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
              time = "visit", treatment = "trt")
easte(fit)
#> Separable treatment effects averaged over covariates (population level)
#> Average over 1000 rows of the covariates of the fitted data
#> Thresholds: y1 <= 0, y2 <= 0; Wald 95% confidence intervals
#> Standard errors include the sampling variability of the covariates
#> 
#> End time K = 10:
#>        Estimate Std. Error Lower 95% Upper 95%
#> ESTE1  -0.06862    0.02125   -0.1103  -0.02696
#> ESTE2  -0.04283    0.01739  -0.07691 -0.008757
#> ETTE   -0.05573    0.01040  -0.07612  -0.03534
#> EASTE1  0.06868    0.02122   0.02710    0.1103
#> EASTE2  0.04283    0.01739  0.008757   0.07691
#> DEASTE  0.02585    0.03276  -0.03835   0.09005
#> H0: DEASTE >= 0 vs H1: DEASTE < 0: z = 0.7892, p-value = 0.785

# averaging over draws from a model for the covariates instead
set.seed(1)
draws <- data.frame(x1 = rnorm(2000), x2 = rbinom(2000, 1, 0.4))
easte(fit, newdata = draws)
#> Separable treatment effects averaged over covariates (population level)
#> Average over 2000 rows of the supplied covariate values
#> Thresholds: y1 <= 0, y2 <= 0; Wald 95% confidence intervals
#> Standard errors exclude the sampling variability of the covariates
#> 
#> End time K = 10:
#>        Estimate Std. Error Lower 95% Upper 95%
#> ESTE1  -0.06764    0.02111   -0.1090  -0.02627
#> ESTE2  -0.04252    0.01726  -0.07635 -0.008702
#> ETTE   -0.05508    0.01033  -0.07533  -0.03484
#> EASTE1  0.06782    0.02099   0.02668    0.1090
#> EASTE2  0.04252    0.01726  0.008702   0.07635
#> DEASTE  0.02530    0.03245  -0.03831   0.08890
#> H0: DEASTE >= 0 vs H1: DEASTE < 0: z = 0.7795, p-value = 0.782
```
