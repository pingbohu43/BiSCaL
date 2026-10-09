# Simulation studies of the paper

This vignette reproduces the simulation studies of Section 5.2 of Hu and
Yi (2026) with `BiSCaL`. It describes the designs, computes the true
values of the effects, analyzes one simulated data set, and summarizes a
full run of the studies (100 replicates of each design).

``` r

library(BiSCaL)
```

## Designs

Each simulated data set has n = 5000 subjects followed at a baseline
visit and K = 20 follow-up visits. There is one standard normal
covariate X, and a treatment assigned with probability 0.5. The binary
responses follow models (12)–(14) of the paper. All coefficients equal
0.5 except the treatment coefficients (\beta\_{1,3,1}, \beta\_{2,3,1}):

- case 1: (5, 0.1), the treatment mainly affects the first response;
- case 2: (0, 0), no treatment effect;
- case 3: (0.1, 5), the treatment mainly affects the second response.

``` r

biscal_coefs(case = 1)$beta1
#>               1
#> (Intercept) 0.5
#> y1[k-1]=1   0.5
#> y2[k-1]=1   0.5
#> trt         5.0
#> x           0.5
```

[`simulate_biscal()`](https://pingbohu43.github.io/BiSCaL/reference/simulate_biscal.md)
generates complete data (Section 5.2.2), or data with censoring and
intermittent missingness (`missing = TRUE`, Section 5.2.3). In the
latter, subject i is followed up to visit C_i \sim \mathrm{Binomial}(20,
\sigma(1.5 + x_i)). At each follow-up visit, the two responses are
missing with probabilities \sigma(x_i) and \sigma(x_i - 0.5).

``` r

d <- simulate_biscal(n = 5000, K = 20, case = 3, missing = TRUE, seed = 1)
head(d)
#>   id visit trt          x y1 y2
#> 1  1     0   0 -0.6264538  0  0
#> 2  1     1   0 -0.6264538  0  0
#> 3  1     2   0 -0.6264538  0  1
#> 4  1     3   0 -0.6264538  1  1
#> 5  1     4   0 -0.6264538  1  1
#> 6  1     5   0 -0.6264538  1  1
mean(attr(d, "censoring"))                     # mean censoring time
#> [1] 15.5218
fu <- d[d$visit > 0, ]
c(mean(tapply(is.na(fu$y1), fu$id, mean)),     # mean missingness rates
  mean(tapply(is.na(fu$y2), fu$id, mean)))
#> [1] 0.5018990 0.3963568
```

## True values

The target is DEASTE_K = E\|STE\_{1,K}(X)\| - E\|STE\_{2,K}(X)\| at K =
20, the expectation being over X \sim N(0, 1).
[`biscal_truth()`](https://pingbohu43.github.io/BiSCaL/reference/biscal_truth.md)
computes it by numerical integration:

``` r

truth <- t(sapply(1:3, function(case) biscal_truth(K = 20, case = case)))
rownames(truth) <- paste("case", 1:3)
signif(truth, 4)
#>             ESTE1     ESTE2      ETTE    EASTE1   EASTE2    DEASTE
#> case 1 -0.0099080 -0.003365 -0.006636 0.0099080 0.003393  0.006515
#> case 2  0.0000000  0.000000  0.000000 0.0000000 0.000000  0.000000
#> case 3 -0.0001733 -0.093830 -0.047000 0.0001755 0.093830 -0.093660
```

These are the values T in Tables 5 and 7 of the paper. Case 1 gives
6.515 \times 10^{-3}, against 6.535 \times 10^{-3} in the paper, whose
value was computed by Monte Carlo integration.

## One replicate

Each replicate fits the models, then estimates DEASTE\_{20} as an
average over 10,000 draws from the covariate distribution. It tests H_0:
DEASTE\_{20} \ge 0 against H_1: DEASTE\_{20} \< 0:

``` r

fit <- biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit", treatment = "trt")
set.seed(2)
draws <- data.frame(x = rnorm(10000))
res <- easte(fit, newdata = draws, K = 20)
res$test
#>    K   estimate  std.error statistic      p.value
#> 1 20 -0.1299243 0.01597084 -8.135095 2.058078e-16
```

The estimate from a single data set carries sampling error. The repeated
replicates below show the behavior of the estimator.

The full study repeats this 100 times for each design. The code is in
`data-raw/simulation-study.R` in the package sources. It uses
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html) and
takes about 2 seconds per replicate and core.

``` r

one_rep <- function(i, case, missing) {
  seed <- 100000 * missing + 1000 * case + i
  d <- simulate_biscal(n = 5000, K = 20, case = case, missing = missing, seed = seed)
  fit <- biscal(cbind(y1, y2) ~ x, data = d, id = "id", time = "visit",
                treatment = "trt")
  set.seed(seed + 500000)
  e <- easte(fit, newdata = data.frame(x = rnorm(10000)), K = 20)
  list(coef = coef(fit), se = sqrt(diag(vcov(fit))), deaste = e$test$estimate,
       se_deaste = e$test$std.error, p = e$test$p.value)
}
reps <- parallel::mclapply(1:100, one_rep, case = 1, missing = FALSE, mc.cores = 4)
```

## Results of the full study

The results of the full run are stored in the package:

``` r

sim <- readRDS(system.file("extdata", "simulation-study.rds", package = "BiSCaL"))
sim$settings
#> $n
#> [1] 5000
#> 
#> $K
#> [1] 20
#> 
#> $n_draws
#> [1] 10000
#> 
#> $R
#> [1] 100
```

### Separable effects (Tables 5 and 7)

For each design, the table gives the true value (T) and the average
estimate (AE). It also gives the average standard error (ASE), the
average p-value of the test (AP) and the coverage (CP, %) of the 95%
confidence interval:

``` r

summ_deaste <- function(r) {
  covered <- abs(r$deaste - r$true_deaste) <= qnorm(0.975) * r$se_deaste
  c(T = r$true_deaste, AE = mean(r$deaste), ASE = mean(r$se_deaste),
    AP = mean(r$p), CP = 100 * mean(covered))
}
tab <- t(sapply(sim$results, summ_deaste))
rownames(tab) <- sprintf("case %d, %s", sapply(sim$results, `[[`, "case"),
                         ifelse(sapply(sim$results, `[[`, "missing"),
                                "missingness and censoring", "complete data"))
signif(tab, 3)
#>                                          T        AE      ASE       AP CP
#> case 1, complete data              0.00652  0.006390 0.001290 1.00e+00 93
#> case 2, complete data              0.00000 -0.000392 0.000541 2.72e-01 96
#> case 3, complete data             -0.09370 -0.094500 0.006240 1.39e-47 98
#> case 1, missingness and censoring  0.00652  0.006100 0.002890 9.42e-01 96
#> case 2, missingness and censoring  0.00000 -0.000782 0.001230 3.02e-01 96
#> case 3, missingness and censoring -0.09370 -0.093700 0.013300 1.15e-10 91
```

As in the paper, the sign of the estimated DEASTE\_{20} matches the
design:

- case 1: positive, and the test does not reject H_0;
- case 2: close to zero;
- case 3: negative, and H_0 is rejected.

The coverage of the confidence intervals is close to the nominal 95%.

### Model parameters (Tables 4 and 6)

The same summaries for the coefficients of the regression models; the
labels follow the notation of the paper:

``` r

labels <- c("beta[1,0,1]", "beta[1,1,1]", "beta[1,2,1]", "beta[1,3,1]", "beta[1,4,1]",
            "beta[2,0,1]", "beta[2,1,1]", "beta[2,2,1]", "beta[2,3,1]", "beta[2,4,1]",
            "gamma[1,0,1]", "gamma[1,2,1]", "gamma[2,0,1]", "gamma[2,1,1]", "gamma[2,2,1]")
summ_coef <- function(r) {
  covered <- abs(sweep(r$coef, 2, r$true_coef)) <= qnorm(0.975) * r$se
  data.frame(T = r$true_coef, AE = colMeans(r$coef), ASE = colMeans(r$se),
             CP = 100 * colMeans(covered), row.names = labels)
}
for (r in sim$results) {
  cat(sprintf("\nCase %d, %s\n", r$case,
              if (r$missing) "missingness and censoring" else "complete data"))
  print(signif(summ_coef(r), 3))
}
#> 
#> Case 1, complete data
#>                T    AE    ASE CP
#> beta[1,0,1]  0.5 0.500 0.0255 96
#> beta[1,1,1]  0.5 0.498 0.0240 95
#> beta[1,2,1]  0.5 0.502 0.0240 95
#> beta[1,3,1]  5.0 5.010 0.1010 96
#> beta[1,4,1]  0.5 0.500 0.0119 91
#> beta[2,0,1]  0.5 0.501 0.0231 97
#> beta[2,1,1]  0.5 0.499 0.0238 97
#> beta[2,2,1]  0.5 0.499 0.0178 94
#> beta[2,3,1]  0.1 0.100 0.0175 95
#> beta[2,4,1]  0.5 0.501 0.0086 95
#> gamma[1,0,1] 0.5 0.499 0.0301 94
#> gamma[1,2,1] 0.5 0.498 0.0315 90
#> gamma[2,0,1] 0.5 0.505 0.0492 95
#> gamma[2,1,1] 0.5 0.498 0.0645 96
#> gamma[2,2,1] 0.5 0.504 0.0340 95
#> 
#> Case 2, complete data
#>                T       AE     ASE CP
#> beta[1,0,1]  0.5 4.97e-01 0.01970 99
#> beta[1,1,1]  0.5 5.03e-01 0.01700 98
#> beta[1,2,1]  0.5 5.01e-01 0.01710 98
#> beta[1,3,1]  0.0 1.49e-05 0.01540 96
#> beta[1,4,1]  0.5 4.98e-01 0.00850 92
#> beta[2,0,1]  0.5 5.00e-01 0.01990 94
#> beta[2,1,1]  0.5 5.00e-01 0.01710 95
#> beta[2,2,1]  0.5 4.99e-01 0.01710 91
#> beta[2,3,1]  0.0 1.45e-03 0.01540 93
#> beta[2,4,1]  0.5 5.01e-01 0.00851 93
#> gamma[1,0,1] 0.5 4.98e-01 0.03010 96
#> gamma[1,2,1] 0.5 5.00e-01 0.03160 99
#> gamma[2,0,1] 0.5 4.94e-01 0.04910 93
#> gamma[2,1,1] 0.5 5.07e-01 0.06440 94
#> gamma[2,2,1] 0.5 5.00e-01 0.03390 96
#> 
#> Case 3, complete data
#>                T     AE     ASE CP
#> beta[1,0,1]  0.5 0.4950 0.02240 97
#> beta[1,1,1]  0.5 0.5020 0.01770 94
#> beta[1,2,1]  0.5 0.5040 0.02280 93
#> beta[1,3,1]  0.1 0.0976 0.01720 93
#> beta[1,4,1]  0.5 0.4980 0.00858 96
#> beta[2,0,1]  0.5 0.4980 0.02580 91
#> beta[2,1,1]  0.5 0.5020 0.02410 94
#> beta[2,2,1]  0.5 0.5000 0.02420 95
#> beta[2,3,1]  5.0 5.0100 0.10200 98
#> beta[2,4,1]  0.5 0.4990 0.01200 95
#> gamma[1,0,1] 0.5 0.5020 0.03010 98
#> gamma[1,2,1] 0.5 0.4990 0.03160 96
#> gamma[2,0,1] 0.5 0.5040 0.04920 95
#> gamma[2,1,1] 0.5 0.4970 0.06440 95
#> gamma[2,2,1] 0.5 0.4960 0.03390 93
#> 
#> Case 1, missingness and censoring
#>                T    AE    ASE CP
#> beta[1,0,1]  0.5 0.499 0.0619 94
#> beta[1,1,1]  0.5 0.497 0.0539 95
#> beta[1,2,1]  0.5 0.501 0.0541 94
#> beta[1,3,1]  5.0 5.020 0.2240 95
#> beta[1,4,1]  0.5 0.502 0.0353 99
#> beta[2,0,1]  0.5 0.503 0.0528 96
#> beta[2,1,1]  0.5 0.497 0.0519 97
#> beta[2,2,1]  0.5 0.496 0.0385 94
#> beta[2,3,1]  0.1 0.104 0.0403 97
#> beta[2,4,1]  0.5 0.502 0.0249 93
#> gamma[1,0,1] 0.5 0.504 0.0301 97
#> gamma[1,2,1] 0.5 0.505 0.0316 97
#> gamma[2,0,1] 0.5 0.493 0.0492 97
#> gamma[2,1,1] 0.5 0.502 0.0644 94
#> gamma[2,2,1] 0.5 0.503 0.0339 98
#> 
#> Case 2, missingness and censoring
#>                T      AE    ASE CP
#> beta[1,0,1]  0.5  0.4980 0.0476 96
#> beta[1,1,1]  0.5  0.4980 0.0382 93
#> beta[1,2,1]  0.5  0.5030 0.0385 98
#> beta[1,3,1]  0.0 -0.0071 0.0365 94
#> beta[1,4,1]  0.5  0.4970 0.0251 93
#> beta[2,0,1]  0.5  0.4980 0.0468 96
#> beta[2,1,1]  0.5  0.5030 0.0373 96
#> beta[2,2,1]  0.5  0.4980 0.0371 93
#> beta[2,3,1]  0.0 -0.0012 0.0350 94
#> beta[2,4,1]  0.5  0.5010 0.0246 92
#> gamma[1,0,1] 0.5  0.5020 0.0301 94
#> gamma[1,2,1] 0.5  0.5050 0.0316 93
#> gamma[2,0,1] 0.5  0.5010 0.0492 97
#> gamma[2,1,1] 0.5  0.5000 0.0645 97
#> gamma[2,2,1] 0.5  0.5020 0.0339 90
#> 
#> Case 3, missingness and censoring
#>                T    AE    ASE CP
#> beta[1,0,1]  0.5 0.492 0.0506 93
#> beta[1,1,1]  0.5 0.500 0.0394 97
#> beta[1,2,1]  0.5 0.502 0.0478 94
#> beta[1,3,1]  0.1 0.101 0.0398 90
#> beta[1,4,1]  0.5 0.498 0.0253 94
#> beta[2,0,1]  0.5 0.506 0.0612 97
#> beta[2,1,1]  0.5 0.498 0.0524 94
#> beta[2,2,1]  0.5 0.498 0.0525 94
#> beta[2,3,1]  5.0 5.040 0.2210 94
#> beta[2,4,1]  0.5 0.503 0.0345 95
#> gamma[1,0,1] 0.5 0.499 0.0301 95
#> gamma[1,2,1] 0.5 0.497 0.0315 92
#> gamma[2,0,1] 0.5 0.508 0.0491 96
#> gamma[2,1,1] 0.5 0.497 0.0644 97
#> gamma[2,2,1] 0.5 0.498 0.0339 95
```

The estimates are close to the true values and the coverage is close to
95%. This reproduces the conclusions of Tables 4 and 6.

The numbers differ from the published tables because the replicates use
different random numbers. For the design with missingness and censoring,
the data generation is also implemented differently: here, no visit
after the censoring time C_i is recorded. With 100 replicates, the Monte
Carlo standard error of a coverage near 95% is about 2.2 percentage
points.

## Reference

Hu, P. and Yi, G. Y. (2026). Inference about separable causal effects
with longitudinal bivariate ordinal responses with missingness and
censoring. *Statistics in Medicine*, 45, e70320.
<https://doi.org/10.1002/sim.70320>
