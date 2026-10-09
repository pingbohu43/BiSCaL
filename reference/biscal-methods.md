# Methods for fitted biscal models

Print, summarize and extract information from a model fitted by
[`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md).

## Usage

``` r
# S3 method for class 'biscal'
coef(object, ...)

# S3 method for class 'biscal'
vcov(object, ...)

# S3 method for class 'biscal'
nobs(object, ...)

# S3 method for class 'biscal'
logLik(object, ...)

# S3 method for class 'biscal'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'biscal'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'biscal'
summary(object, ...)

# S3 method for class 'summary.biscal'
print(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
)
```

## Arguments

- object, x:

  A fitted model returned by
  [`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md)
  (or its summary).

- ...:

  Further arguments (unused).

- parm, level:

  Coefficients (names or indices) and confidence level for
  [`confint()`](https://rdrr.io/r/stats/confint.html).

- digits:

  Number of significant digits to print.

- signif.stars:

  Show significance stars in the coefficient tables?

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) returns the stacked
coefficient vector \\\theta = (\beta_1^\top, \beta_2^\top,
\gamma_1^\top, \gamma_2^\top)^\top\\, with names `model:term` (for
responses with more than two categories `model[category]:term`);
[`vcov()`](https://rdrr.io/r/stats/vcov.html) its estimated covariance
matrix; [`confint()`](https://rdrr.io/r/stats/confint.html) Wald
confidence intervals; [`logLik()`](https://rdrr.io/r/stats/logLik.html)
the maximized log-likelihood;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of subjects.
[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.biscal"` whose `coefficients` element is a list of
coefficient tables (estimate, standard error, z value and two-sided
p-value, as in Table 1 of the paper), one per model.

## Examples

``` r
fit <- biscal(cbind(y1, y2) ~ x1 + x2, data = simtrial, id = "id",
              time = "visit", treatment = "trt")
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
head(coef(fit))
#> beta1:(Intercept)   beta1:y1[k-1]=1   beta1:y2[k-1]=1         beta1:trt 
#>        -0.5761739         1.5983598         1.1998988         1.6198973 
#>          beta1:x1          beta1:x2 
#>         0.4016786        -0.5472111 
confint(fit, parm = "beta1:trt")
#>              2.5 %   97.5 %
#> beta1:trt 1.412751 1.827044
```
