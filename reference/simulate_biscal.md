# Simulate longitudinal bivariate responses

Generates data from the models of the paper: baseline responses from
models (13) and (14), then, at each visit \\k = 1, \dots, K\\, the first
response given the previous responses, the treatment and the covariates,
and the second response given the current first response, its own
previous value, the treatment and the covariates (model (12)). By
default this is the design of the simulation studies of the paper
(Section 5.2): one standard normal covariate, a treatment assigned with
probability 0.5 independently of it, and the coefficients of
[`biscal_coefs()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_coefs.md).

## Usage

``` r
simulate_biscal(
  n = 5000,
  K = 20,
  case = 1,
  missing = FALSE,
  coefs = NULL,
  x = NULL,
  seed = NULL
)
```

## Arguments

- n:

  Number of subjects.

- K:

  Number of follow-up visits; visits are numbered `0, 1, ..., K`.

- case:

  Simulation case of the paper (`1`, `2` or `3`), used when `coefs` is
  `NULL`.

- missing:

  Generate censoring and intermittent missingness as in the paper (see
  below)?

- coefs:

  Optional list of true coefficient matrices with the layout of
  [`biscal_coefs()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_coefs.md);
  more columns give responses with more categories, and more covariate
  rows more covariates.

- x:

  Optional covariate values: a numeric vector or an `n`-row matrix with
  one column per covariate of `coefs`. By default independent standard
  normal variables.

- seed:

  Optional random seed. When given, the result is reproducible and the
  random number stream of the session is left unchanged.

## Value

A data frame in long format, sorted by subject and visit, with columns
`id`, `visit`, `trt` (0/1), the covariate(s) (`x` for a single
covariate, otherwise `x1`, `x2`, ... or the column names of `x`), `y1`
and `y2` (integer categories `0, 1, ...`; `NA` when missing). The
attributes `coefs` and `censoring` hold the true coefficients and, with
`missing = TRUE`, the last visit \\C_i\\ of each subject.

## Missingness and censoring

With `missing = TRUE` the data follow the design of Section 5.2.3 of the
paper. Subject \\i\\ is followed up to visit \\C_i \sim
\mathrm{Binomial}(K, \sigma(1.5 + x_i))\\ (no visit after \\C_i\\ is
recorded), and at each visit \\1 \le k \le C_i\\ the first and second
responses are missing independently with probabilities \\\sigma(x_i)\\
and \\\sigma(x_i - 0.5)\\, where \\\sigma\\ is the logistic function and
\\x_i\\ the first covariate. Baseline responses are always observed.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

[`biscal_coefs()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_coefs.md),
[`biscal_truth()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_truth.md),
[`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md).

## Examples

``` r
d <- simulate_biscal(n = 300, K = 5, case = 1, missing = TRUE, seed = 1)
head(d)
#>   id visit trt          x y1 y2
#> 1  1     0   1 -0.6264538  1  1
#> 2  1     1   1 -0.6264538  1 NA
#> 3  1     2   1 -0.6264538  1  0
#> 4  2     0   1  0.1836433  1  1
#> 5  2     1   1  0.1836433 NA  1
#> 6  2     2   1  0.1836433  1 NA
fit <- biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit",
              treatment = "trt")
coef(fit)
#>  beta1:(Intercept)    beta1:y1[k-1]=1    beta1:y2[k-1]=1          beta1:trt 
#>         0.49007194         0.22810723         0.53198560         4.32442076 
#>            beta1:x  beta2:(Intercept)      beta2:y1[k]=1    beta2:y2[k-1]=1 
#>         0.44911363         0.70403351         0.70138188         0.05337212 
#>          beta2:trt            beta2:x gamma1:(Intercept)           gamma1:x 
#>         0.14078057         0.59522238         0.35802987         0.61659751 
#> gamma2:(Intercept)     gamma2:y1[0]=1           gamma2:x 
#>         0.48692605         0.52188554         0.31281004 
```
