# True effects of a simulation design

Computes the true separable and total treatment effects of a simulation
design with a single standard normal covariate (such as the designs of
the paper), either at given covariate values or averaged over the
covariate distribution, the targets of
[`aste()`](https://pingbohu43.github.io/BiSCaL/reference/aste.md) and
[`easte()`](https://pingbohu43.github.io/BiSCaL/reference/easte.md).
Population-level values are computed by numerical integration.

## Usage

``` r
biscal_truth(K = 20, case = 1, coefs = NULL, y1 = 0, y2 = 0, x = NULL)
```

## Arguments

- K:

  End time point.

- case, coefs:

  The design, as in
  [`simulate_biscal()`](https://pingbohu43.github.io/BiSCaL/reference/simulate_biscal.md).

- y1, y2:

  Thresholds of the generalized odds ratio, as codes `0, 1, ...` of the
  response categories (not the highest).

- x:

  Optional numeric vector of covariate values. If `NULL` (default),
  population-level values under \\X \sim N(0, 1)\\ are returned.

## Value

With `x = NULL`, a named vector with the true `ESTE1`, `ESTE2`, `ETTE`,
`EASTE1`, `EASTE2` and `DEASTE`. Otherwise a data frame with the true
`STE1`, `STE2`, `TTE`, `ASTE1`, `ASTE2` and `DASTE` at each value of
`x`.

## References

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, **45**, e70320.
[doi:10.1002/sim.70320](https://doi.org/10.1002/sim.70320)

## See also

[`simulate_biscal()`](https://pingbohu43.github.io/BiSCaL/reference/simulate_biscal.md),
[`biscal_coefs()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_coefs.md).

## Examples

``` r
# true DEASTE of the three cases (Table 5 of the paper)
sapply(1:3, function(case) biscal_truth(K = 20, case = case)[["DEASTE"]])
#> [1]  0.006515171  0.000000000 -0.093659166
biscal_truth(K = 20, case = 1, x = c(-1, 0, 1))
#>    x         STE1         STE2          TTE       ASTE1       ASTE2       DASTE
#> 1 -1 -0.013793451 -0.002554046 -0.008173749 0.013793451 0.002554046 0.011239405
#> 2  0 -0.010164266 -0.003870253 -0.007017260 0.010164266 0.003870253 0.006294012
#> 3  1 -0.006131398 -0.004093514 -0.005112456 0.006131398 0.004093514 0.002037885
```
