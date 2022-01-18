# Example 1

This tutorial is devoted to comparing the computational performance of
the different corrections on a set of artificial data. We generate
simulated data from a logistic regression model, with correlated
covariates and regression coefficient sampled from a standard gaussian
distribution.

``` r
set.seed(1991)
n <- 50
p <- 10

SS <- matrix(0.2, p - 1, p - 1)
diag(SS) <- 1

X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1) %*% chol(SS))
betas <- runif(p)
y <- rbinom(n, 1, prob = plogis(X %*% betas))
```

We estimate different corrections and compare their empirical
performance in terms of timing. Recalling that penalized likelihood
optimization under the Diaconis-Ylvisaker conjugate prior induces a
binomial likelihood with pseudo-counts, we can rely on different
available optimization methods for logistic regression models. We found
that quasi-Newton methods have better numerical performance; refer to
Nocedal and Wright (2006) for additional details. A practical
implementation for logistic regression optimization via L-BFGS is
provided in the function `fastLR` available in the package
`RcppNumerical`

``` r
library(brglm2)
m <- nrow(X)
y_dy <- p / (p + m) * 0.5 + m / (p + m) * y

# DY ESTIMATE
t0 <- Sys.time()
fit_dy <- glm(y_dy ~ X - 1, family = binomial("logit"))
t1 <- Sys.time()
elapsed_dy <- t1 - t0

library(RcppNumerical)
t0 <- Sys.time()
fit_fast_dy <- fastLR(X, y_dy)
t1 <- Sys.time()
elapsed_fast_dy <- t1 - t0

# FIRTH (1993)
t0 <- Sys.time()
fit_firth <- glm(y ~ -1 + X,
  family = binomial("logit"),
  method = "brglmFit", type = "AS_mean"
)
t1 <- Sys.time()
elapsed_firth <- t1 - t0

# KENNE PAGUI ET AL. (2017)
t0 <- Sys.time()
fit_kp <- glm(y ~ -1 + X,
  family = binomial("logit"),
  method = "brglmFit", type = "AS_median"
)
t1 <- Sys.time()
elapsed_kp <- t1 - t0
```

Timing comparison are provided in the following table.

|                                         | x          |
|:----------------------------------------|:-----------|
| DY                                      | 0.004 secs |
| DY - fast implementation                | 0.002 secs |
| Firth (1993)                            | 0.010 secs |
| Kenne Pagui, Salvan, and Sartori (2017) | 0.060 secs |

Computational costs notably increases with *n* and *p*, keeping the same
ratio *p*/*n* = 0.2. For examples, with *n* = 500 and *p* = 100

|                                         | x           |
|:----------------------------------------|:------------|
| DY                                      | 0.039 secs  |
| DY - fast implementation                | 0.003 secs  |
| Firth (1993)                            | 1.374 secs  |
| Kenne Pagui, Salvan, and Sartori (2017) | 33.086 secs |

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Firth1993" class="csl-entry">

Firth, David. 1993. “<span class="nocase">Bias reduction of maximum
likelihood estimates</span>.” *Biometrika* 80 (1): 27–38.

</div>

<div id="ref-Pagui2017" class="csl-entry">

Kenne Pagui, E. C., A. Salvan, and N. Sartori. 2017. “<span
class="nocase">Median bias reduction of maximum likelihood
estimates</span>.” *Biometrika* 104 (4): 923–38.

</div>

<div id="ref-Nocedal2006" class="csl-entry">

Nocedal, Jorge, and Stephen Wright. 2006. *Conjugate Gradient Methods*.
New York, NY: Springer New York.

</div>

</div>
