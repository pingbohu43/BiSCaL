# Creates data/simtrial.rda, the simulated example data set (run from the
# package root with the package loaded, e.g. after devtools::load_all()).
# Binary responses and two baseline covariates; the treatment mainly acts on
# the first response. Subjects are censored as in Section 5.2.3 of the paper;
# intermittent missingness is lighter than in the paper's simulation design.
mk <- function(v, terms) matrix(v, ncol = 1L, dimnames = list(terms, "1"))
coefs <- list(
  beta1 = mk(c(-0.5, 1.5, 1.0, 1.5, 0.4, -0.3),
             c("(Intercept)", "y1[k-1]=1", "y2[k-1]=1", "trt", "x1", "x2")),
  beta2 = mk(c(-0.5, 1.0, 1.5, 0.2, 0.3, 0.2),
             c("(Intercept)", "y1[k]=1", "y2[k-1]=1", "trt", "x1", "x2")),
  gamma1 = mk(c(0, 0.5, -0.3), c("(Intercept)", "x1", "x2")),
  gamma2 = mk(c(0, 1.0, 0.4, 0.2), c("(Intercept)", "y1[0]=1", "x1", "x2"))
)
n <- 1000
K <- 10
set.seed(2026)
X <- cbind(x1 = round(rnorm(n), 3), x2 = rbinom(n, 1, 0.4))
full <- simulate_biscal(n = n, K = K, coefs = coefs, x = X, seed = 2027)

# censoring: last visit C ~ Binomial(K, plogis(1.5 + x1)); intermittent
# missingness at follow-up visits with probabilities plogis(-1.5 + 0.5 x1)
# (first response) and plogis(-2 + 0.5 x1) (second response)
set.seed(2028)
cens <- rbinom(n, K, plogis(1.5 + X[, "x1"]))
x1 <- full$x1
fu <- full$visit > 0
miss1 <- fu & runif(nrow(full)) < plogis(-1.5 + 0.5 * x1)
miss2 <- fu & runif(nrow(full)) < plogis(-2 + 0.5 * x1)
full$y1[miss1] <- NA
full$y2[miss2] <- NA
simtrial <- full[full$visit <= cens[full$id], ]
rownames(simtrial) <- NULL
attr(simtrial, "censoring") <- NULL
attr(simtrial, "coefs") <- coefs
save(simtrial, file = "data/simtrial.rda", compress = "xz")
