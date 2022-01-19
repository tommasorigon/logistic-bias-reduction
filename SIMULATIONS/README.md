# Simulation studies

In this tutorial we focus on different simulation studies to evaluate
the performance of the proposed DY correction against competitors.

## First example

The first example is taken from Appendix D Sur and Candès (2019) and the
Supplementary Materialds of Kosmidis and Firth (2021), and focuses on a
setting with *n* = 1000 observations from a logisitc regression model
with *p* = 200, 25 of coefficients equal to 10, 25 of coefficients equal
to  − 10 and the remaining set as 0.

``` r
set.seed(1)
n <- 2000
p <- 400
beta <- c(rep(10, p / 8), rep(-10, p / 8), rep(0, 3 * p / 4))
X <- matrix(rnorm(n * p, 0, sqrt(1 / n)), n, p)
y <- rbinom(n, 1, plogis(X %*% beta))
fit_mle <- glm(y ~ X - 1, family = binomial("logit"))
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




library(RcppNumerical)
t0 <- Sys.time()
fit_fast_dy <- fastLR(X, y_dy)
t1 <- Sys.time()
elapsed_fast_dy <- t1 - t0
y_clogg <- p / (p + m) * mean(y) + m / (p + m) * y

# DY ESTIMATE
fit_fast_clogg <- fastLR(X, y_clogg)

# FIRTH (1993)
t0 <- Sys.time()
fit_firth <- glm(y ~ -1 + X,
  family = binomial("logit"),
  method = "brglmFit", type = "AS_mean"
)
t1 <- Sys.time()
elapsed_firth <- t1 - t0
```

<img src="figs/coef.png" style="display: block; margin: auto;" />

We evaluate bias and rmse across 1000 replications of this scenario.
Computations takes roughly 5 hours on a 2020 Macbook Pro with M1
processor (`aarch64-apple-darwin20`) running R 4.1.1 linked with
`openblas`. Results are stored in the file `sur-candes.RData` and can be
reproduced running the script
[`sur-candes.R`](https://raw.githubusercontent.com/tommasorigon/logistic-bias-reduction/main/SIMULATIONS/sur-candes.R)

<img src="figs/boxpl-1.png" style="display: block; margin: auto;" />

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Kosmidis2021" class="csl-entry">

Kosmidis, Ioannis, and David Firth. 2021. “<span
class="nocase">Jeffreys-prior penalty, finiteness and shrinkage in
binomial-response generalized linear models</span>.” *Biometrika* 108
(1): 71–82. <https://doi.org/10.1093/biomet/asaa052>.

</div>

<div id="ref-Nocedal2006" class="csl-entry">

Nocedal, Jorge, and Stephen Wright. 2006. *Conjugate Gradient Methods*.
New York, NY: Springer New York.

</div>

<div id="ref-Sur2019" class="csl-entry">

Sur, P, and E. J. Candès. 2019. “A Modern Maximum-Likelihood Theory for
High-Dimensional Logistic Regression.” *Proc. Nat. Acad. Sci.* 116:
14516–25.

</div>

</div>
