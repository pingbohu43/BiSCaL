# Coefficients of the simulation designs of the paper

Returns the true coefficients of the three cases of the simulation
studies in Section 5.2 of the paper: binary responses, one standard
normal covariate, all coefficients equal to 0.5 except the treatment
coefficients \\(\beta\_{1,3,1}, \beta\_{2,3,1})\\, which are \\(5,
0.1)\\ in case 1 (the treatment mainly affects the first response),
\\(0, 0)\\ in case 2 (no treatment effect) and \\(0.1, 5)\\ in case 3
(the treatment mainly affects the second response).

## Usage

``` r
biscal_coefs(case = 1)
```

## Arguments

- case:

  `1`, `2` or `3`.

## Value

A list of one-column matrices `beta1`, `beta2`, `gamma1` and `gamma2`,
laid out as the `coefficients` of a
[`biscal()`](https://pingbohu43.github.io/BiSCaL/reference/biscal.md)
fit: terms in rows (intercept, the dummy-coded previous responses, the
treatment `trt` and the covariate `x`, in the order of models (12)-(14)
of the paper) and outcome categories in columns.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

[`simulate_biscal()`](https://pingbohu43.github.io/BiSCaL/reference/simulate_biscal.md),
[`biscal_truth()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_truth.md).

## Examples

``` r
biscal_coefs(1)
#> $beta1
#>               1
#> (Intercept) 0.5
#> y1[k-1]=1   0.5
#> y2[k-1]=1   0.5
#> trt         5.0
#> x           0.5
#> 
#> $beta2
#>               1
#> (Intercept) 0.5
#> y1[k]=1     0.5
#> y2[k-1]=1   0.5
#> trt         0.1
#> x           0.5
#> 
#> $gamma1
#>               1
#> (Intercept) 0.5
#> x           0.5
#> 
#> $gamma2
#>               1
#> (Intercept) 0.5
#> y1[0]=1     0.5
#> x           0.5
#> 
```
