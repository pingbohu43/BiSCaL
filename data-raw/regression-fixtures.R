# Creates tests/testthat/fixtures/original-code.rds: reference values computed
# with the authors' original research code (the simulation functions in
# simulation/FunctionSets.R and the effect functions of the data analysis in
# real data analysis/EssentialFunctions.R). Run from the package root:
#   BISCAL_CODE_DIR=/path/to/code Rscript data-raw/regression-fixtures.R
# The original functions are used unchanged, except that the sample sizes
# hard-coded in SimulationFunction() (5000 subjects, 10000 Monte Carlo draws)
# are made arguments so that the reference run takes minutes rather than days.
code_dir <- Sys.getenv("BISCAL_CODE_DIR")
stopifnot(nzchar(code_dir), dir.exists(code_dir))
suppressPackageStartupMessages({library(purrr); library(boot); library(pracma)})

# --- simulation functions (LaplacesDemon::invlogit defined with its formula) ---
invlogit <- function(x) { InvLogit <- 1 / {1 + exp(-x)}; return(InvLogit) }
src <- readLines(file.path(code_dir, "simulation", "FunctionSets.R"))
src <- src[!grepl("^library\\(", src)]
orig <- new.env()
assign("invlogit", invlogit, envir = orig)
eval(parse(text = src), envir = orig)

fn <- deparse(orig$SimulationFunction)
fn <- gsub("5000", "nSubj", fn, fixed = TRUE)
fn <- gsub("10000", "nDraw", fn, fixed = TRUE)
fn[1] <- sub("function (i, alpha, beta, K) ", "function (i, alpha, beta, K, nSubj, nDraw) ",
             fn[1], fixed = TRUE)
at <- grep("result = list()", fn, fixed = TRUE)
fn <- append(fn, paste("assign('simdata', list(x = preteatmentCovariate, trt = treatment,",
                       "Y1 = Y1k, Y2 = Y2k, draws = preteatmentCovariateAlternative),",
                       "envir = globalenv())"), after = at - 1)
SimSmall <- eval(parse(text = fn), envir = orig)
environment(SimSmall) <- orig

K <- 5; nS <- 400; nD <- 300
res <- suppressWarnings(SimSmall(1, c(0.5, 0.5, 0.5, 5, 0.5), c(0.5, 0.5, 0.5, 0.1, 0.5),
                                 K, nS, nD))
long <- data.frame(id = rep(seq_len(nS), each = K + 1), visit = rep(0:K, nS),
                   trt = rep(simdata$trt, each = K + 1), x = rep(simdata$x, each = K + 1),
                   y1 = as.vector(t(simdata$Y1)), y2 = as.vector(t(simdata$Y2)))
sim_complete <- list(long = long, K = K, x = simdata$x, draws = simdata$draws,
                     coef = res[[1]], deaste = res[[2]], ci = res[[3]],
                     var = as.numeric(res[[4]]), vcov = res[[5]])

# true DEASTE at K = 20 for cases 1 and 3, by numerical integration over N(0, 1)
truth <- sapply(list(case1 = c(5, 0.1), case3 = c(0.1, 5)), function(tc) {
  th <- c(0.5, 0.5, rep(0.5, 3), c(0.5, 0.5, 0.5, tc[1], 0.5), c(0.5, 0.5, 0.5, tc[2], 0.5))
  h <- function(x) sapply(x, function(xx) -orig$SimulationSTE2MinusSTE1(20, xx, th))
  integrate(function(x) h(x) * dnorm(x), -10, 10, rel.tol = 1e-10, subdivisions = 1000L)$value
})

# --- effect functions of the data analysis (K = 64, 8 covariates) ---------------
inv.logit <- boot::inv.logit
txt <- readLines(file.path(code_dir, "real data analysis", "EssentialFunctions.R"),
                 warn = FALSE)
ef <- new.env(); assign("inv.logit", inv.logit, envir = ef)
eval(parse(text = txt), envir = ef)
# variant in which the model for Y_{1,1} uses y_{2,0}, as model (12) of the
# paper (the analysis code sets y_{2,0} to 0 in that model)
txt2 <- sub("OutputP11IncludingParameter(y11, i, t11, x1, theta)",
            "OutputP1kIncludingParameter(y11, i, j, t11, x1, theta)", txt, fixed = TRUE)
ef2 <- new.env(); assign("inv.logit", inv.logit, envir = ef2)
eval(parse(text = txt2), envir = ef2)

# Table 1 estimates, in the order of the authors' 43-dimensional theta
b1 <- c(-3.259, 0.178, -0.448, 0.112, 0.444, 3.556e-4, 0.104, 0.236, -0.545, 0.024, 3.443, 0.102)
b2 <- c(-6.137, 2.664, 0.064, -0.015, -0.077, -0.139, -0.012, -0.112, -0.054, 0.027, 3.021, 0.057)
g1 <- c(-0.654, -1.027, 0.137, 0.552, -0.317, 0.258, 0.145, -1.224, 0.030)
g2 <- c(-17.889, -1.256, -1.500, -0.259, 0.983, 0.671, -0.152, 14.787, -0.013, 2.385)
theta <- c(b1, b2, g1, g2)
xmin <- c(1.6, 0, 0, 0, 0, 0, 0, 35)
xmax <- c(4.5, 1, 1, 1, 1, 1, 1, 60)
xn <- c("BMI", "EDUC", "EVDNK", "EVSMK", "PHYSA", "PHYSB", "RACE", "V5AGE")
mk <- function(v, terms) matrix(v, ncol = 1, dimnames = list(terms, "1"))
coefs <- list(beta1 = mk(b1[c(1, 11, 12, 2, 3:10)], c("(Intercept)", "hdl[k-1]=1", "ldl[k-1]=1", "trt", xn)),
              beta2 = mk(b2[c(1, 12, 11, 2, 3:10)], c("(Intercept)", "hdl[k]=1", "ldl[k-1]=1", "trt", xn)),
              gamma1 = mk(g1, c("(Intercept)", xn)),
              gamma2 = mk(g2[c(1, 10, 2:9)], c("(Intercept)", "hdl[0]=1", xn)))
ste <- function(e) rbind(min = c(e$STE1IncludingParameter(64, xmin, theta), e$STE2IncludingParameter(64, xmin, theta)),
                         max = c(e$STE1IncludingParameter(64, xmax, theta), e$STE2IncludingParameter(64, xmax, theta)))
realdata <- list(coefs = coefs, xmin = xmin, xmax = xmax, ste_paper_model = ste(ef2),
                 ste_as_run = ste(ef))

dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
saveRDS(list(sim_complete = sim_complete, truth = truth, realdata = realdata),
        "tests/testthat/fixtures/original-code.rds", version = 2)
