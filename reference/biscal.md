# Fit the regression models for separable treatment effects

Fits, by maximum likelihood, the four baseline-category logistic
regression models through which the separable treatment effects of Hu
and Yi (2026) are identified: the transition models for the two
responses (model (12) of the paper) and the models for their baseline
values (models (13) and (14)). Intermittent missingness and
right-censoring are handled as in Section 4 of the paper.

## Usage

``` r
biscal(formula, data, id, time, treatment, control = list())
```

## Arguments

- formula:

  A formula `cbind(y1, y2) ~ x1 + x2 + ...`. The left-hand side names
  the two responses, the one measured first at each visit first. They
  may be numeric (the sorted distinct values are the ordered
  categories), logical, or factors (ordered or not; the level order is
  used). The right-hand side lists baseline covariates, which may
  include factors, interactions and transformations such as
  [`scale()`](https://rdrr.io/r/base/scale.html) or
  [`poly()`](https://rdrr.io/r/stats/poly.html) (evaluated on new
  covariate values with the parameters of the fitted data); use `~ 1`
  for none. [`offset()`](https://rdrr.io/r/stats/offset.html) terms are
  not supported.

- data:

  A data frame in long format: one row per subject and visit.

- id, time, treatment:

  Names (strings) of the columns of `data` holding the subject
  identifier, the visit index (`0` for the baseline visit, `1, ..., K`
  for the scheduled follow-up visits) and the treatment. The treatment
  may be coded 0/1, logical, or a two-level factor whose first level is
  the control group. The treatment and covariates must be constant
  within subjects.

- control:

  A list of convergence settings: `epsilon` (relative change in the
  log-likelihood, default `1e-8`) and `maxit` (maximum number of
  iterations, default `25`; at least 100 are allowed for responses with
  more than two categories).

## Value

An object of class `"biscal"`, a list with components including

- `coefficients`:

  List of the four coefficient matrices `beta1`, `beta2`, `gamma1` and
  `gamma2` (terms in rows, categories \\d = 1, \dots, L - 1\\ in
  columns).

- `vcov`:

  Block-diagonal covariance matrix of all coefficients (inverse observed
  information), in the order of
  [`coef.biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md).

- `K`:

  The last visit index \\K\\, the default end time point of the effects.

- `n`, `nobs_models`, `loglik`:

  Number of subjects, number of observations used by each model and the
  maximized log-likelihood.

- `X`, `A`, `covariate_data`:

  Subject-level covariate matrix (used by
  [`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md) by
  default), treatment and covariate values.

and information on the data used by the methods.

## Model

For subject \\i\\ at visit \\k = 0, 1, \dots, K\\, let \\Y\_{1,k}\\ and
\\Y\_{2,k}\\ be the two responses, with \\Y\_{1,k}\\ measured before
\\Y\_{2,k}\\; let \\A\\ be the treatment and \\X\\ the baseline
covariates. With \\D_j(\cdot)\\ the dummy (indicator) coding of a
previous response, the models are, for every category \\d \ge 1\\:

- `beta1` (model (12), \\j = 1\\): \\\log\\P(Y\_{1,k} = d \mid
  \cdot)/P(Y\_{1,k} = 0 \mid \cdot)\\ = \beta\_{1,0,d} +
  \beta\_{1,1,d}^\top D_1(y\_{1,k-1}) + \beta\_{1,2,d}^\top
  D_2(y\_{2,k-1}) + \beta\_{1,3,d}\\ a + \beta\_{1,4,d}^\top x\\;

- `beta2` (model (12), \\j = 2\\): the same with \\D_1(y\_{1,k})\\, the
  first response at the same visit, in place of \\D_1(y\_{1,k-1})\\;

- `gamma1` (model (13)): \\Y\_{1,0}\\ given \\x\\;

- `gamma2` (model (14)): \\Y\_{2,0}\\ given \\D_1(y\_{1,0})\\ and \\x\\.

The transition models are common to all visits \\k \ge 1\\. Binary
responses give ordinary logistic regressions. Every category of each
response must occur as an outcome of each model (for example, a response
cannot be constant at baseline), and the visits must be numbered
consecutively so that transitions between visits \\k - 1\\ and \\k\\ are
observed.

## Missing data and censoring

Responses that are not observed are coded `NA`, or their rows are
omitted (for example, all visits after a subject is censored). Following
Section 4.2 of the paper, the transition of the first response at visit
\\k\\ contributes to the likelihood when \\Y\_{1,k}\\, \\Y\_{1,k-1}\\
and \\Y\_{2,k-1}\\ are observed, and that of the second response when
\\Y\_{2,k}\\, \\Y\_{2,k-1}\\ and \\Y\_{1,k}\\ are observed (the
indicators \\m'\_{i,j,k}\\ of the paper). This is valid when missingness
and censoring satisfy Conditions AC1 and AC2 of the paper (for example,
when they depend only on the treatment, the covariates and earlier
responses). The baseline models use the subjects whose baseline
responses are observed.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

[biscal-methods](https://pingbohu43.github.io/BiSCaL/reference/biscal-methods.md)
for [`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html) and other methods;
[`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) and
[`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md) for
the separable treatment effects and tests.

## Examples

``` r
fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
              time = "visit", treatment = "trt")
fit
#> Separable treatment effects for longitudinal bivariate responses (BiSCaL)
#> 
#> Call:
#> biscal(formula = cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id", 
#>     time = "visit", treatment = "trt")
#> 
#> Subjects: 1000 (treated 500, control 500); visits 0 to 10 (K = 10)
#> Treatment: trt (control = 0, treated = 1)
#> Responses: y1 [0 < 1] then y2 [0 < 1]
#> Follow-up: mean last visit 7.77; missing after baseline: y1 20.9%, y2 14.2%
#> Observations used: transitions y1 4511, y2 4700; baseline y1 1000, y2 1000
#> Log-likelihood: -4792.67 (df = 19)
#> 
#> Coefficients:
#> Transition model for y1 (model (12), j = 1):
#> (Intercept)   y1[k-1]=1   y2[k-1]=1         trt          x1          x2 
#>     -0.5762      1.5984      1.1999      1.6199      0.4017     -0.5472 
#> Transition model for y2 (model (12), j = 2):
#> (Intercept)     y1[k]=1   y2[k-1]=1         trt          x1          x2 
#>    -0.45450     0.88056     1.57409     0.22577     0.26392     0.08216 
#> Baseline model for y1 (model (13)):
#> (Intercept)          x1          x2 
#>    -0.05444     0.47013    -0.24870 
#> Baseline model for y2 (model (14)):
#> (Intercept)     y1[0]=1          x1          x2 
#>     -0.0786      0.9546      0.4048      0.3028 
summary(fit)
#> Separable treatment effects for longitudinal bivariate responses (BiSCaL)
#> 
#> Call:
#> biscal(formula = cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id", 
#>     time = "visit", treatment = "trt")
#> 
#> Subjects: 1000 (treated 500, control 500); visits 0 to 10 (K = 10)
#> Treatment: trt (control = 0, treated = 1)
#> Responses: y1 [0 < 1] then y2 [0 < 1]
#> Follow-up: mean last visit 7.77; missing after baseline: y1 20.9%, y2 14.2%
#> Observations used: transitions y1 4511, y2 4700; baseline y1 1000, y2 1000
#> Log-likelihood: -4792.67 (df = 19)
#> 
#> Transition model for y1 (model (12), j = 1):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.57617    0.10007  -5.758 8.54e-09 ***
#> y1[k-1]=1    1.59836    0.09456  16.903  < 2e-16 ***
#> y2[k-1]=1    1.19990    0.09775  12.275  < 2e-16 ***
#> trt          1.61990    0.10569  15.327  < 2e-16 ***
#> x1           0.40168    0.05598   7.176 7.19e-13 ***
#> x2          -0.54721    0.09368  -5.842 5.17e-09 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Transition model for y2 (model (12), j = 2):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.45450    0.09783  -4.646 3.39e-06 ***
#> y1[k]=1      0.88056    0.09867   8.924  < 2e-16 ***
#> y2[k-1]=1    1.57409    0.08364  18.821  < 2e-16 ***
#> trt          0.22577    0.08574   2.633  0.00846 ** 
#> x1           0.26392    0.04791   5.509 3.61e-08 ***
#> x2           0.08216    0.08254   0.995  0.31953    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Baseline model for y1 (model (13)):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.05444    0.08375  -0.650   0.5157    
#> x1           0.47013    0.06911   6.803 1.03e-11 ***
#> x2          -0.24870    0.13322  -1.867   0.0619 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Baseline model for y2 (model (14)):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -0.07860    0.10735  -0.732   0.4640    
#> y1[0]=1      0.95459    0.14091   6.774 1.25e-11 ***
#> x1           0.40480    0.07319   5.531 3.19e-08 ***
#> x2           0.30276    0.14049   2.155   0.0312 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
