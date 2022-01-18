# Birthweight Study

*(This tutorial illustrates main results of the analysis in an
easy-to-read format. Complete code associated with this page is
available in the file*
[birthweight.Rmd](https://raw.githubusercontent.com/tommasorigon/logistic-bias-reduction/main/BIRTHWEIGHT/birthweight.Rmd)),

This tutorial is devoted to replicating simulations from Kosmidis, Kenne
Pagui, and Sartori (2020), available in Table 2 of the article Data
comprises *n* = 100 births and the binary outcome of interest is a
dichotomization of infant birthweight (low against normal).

``` r
library(MASS)
## PREPARE THE DATASET
bwt <- with(birthwt, {
  age <- age
  racewhite <- ifelse(race == 1, 1, 0)
  smoke <- smoke
  ptl <- ifelse(ptl > 0, 1, 0)
  ptd <- factor(ptl > 0)
  ht <- ht
  loglwt <- log(lwt)
  data.frame(normwt = 1 - low, age, racewhite, smoke, ptl, ht, loglwt, ftv)
})
bwt <- subset(bwt, subset = (ftv == 0), select = -c(ftv))
```

A logistic regression is employed with covariates pertaining the mother,
namely: age, race, smoking habitudes, history of premature labour and
hypertension and the logarithm of the weight at birth of the mother.

Similarly to the `endometrial` [case study](../ENDOMETRIAL/), we compare
the proposed Diaconis-Ylvisaker regularization with Clogg et al. (1991),
Firth (1993) and Kenne Pagui, Salvan, and Sartori (2017), relying on the
package `brglm2` (Kosmidis and Firth 2021).

``` r
## ESTIMATES
X <- model.matrix(normwt ~ age + racewhite + smoke + ptl + ht + loglwt, data = bwt)
p <- ncol(X)
n <- m <- nrow(X)

# MLE
ml_fit <- glm(normwt ~ age + racewhite + smoke + ptl + ht + loglwt, family = binomial, data = bwt)

# DY
y_dy <- p / (p + m) * 0.5 + m / (p + m) * bwt$normwt
dy_fit <- glm(y_dy ~ age + racewhite + smoke + ptl + ht + loglwt, family = binomial, data = bwt)

# CLOGG ET AL. (1991)
ybar <- mean(bwt$normwt)
y_clogg <- p / (p + m) * ybar + m / (p + m) * bwt$normwt
clogg_fit <- glm(y_clogg ~ age + racewhite + smoke + ptl + ht + loglwt, family = binomial, data = bwt)

# FIRTH (1993)
br_fit <- update(ml_fit, method = "brglmFit", type = "AS_mean", data = bwt)

# KENNE PAGUI ET AL. (2017)
mbr_fit <- update(ml_fit, method = "brglmFit", type = "AS_median", data = bwt)
```

Estimated regression coefficients are reported in the following table.
Note that maximum-likelihood estimate exists finite in this example.

|                                         |              |              |             |              |              |              |             |
|:----------------------------------------|:-------------|:-------------|:------------|:-------------|:-------------|:-------------|:------------|
| MLE                                     | -8.5 (5.83)  | -0.07 (0.05) | 0.69 (0.57) | -0.56 (0.58) | -1.6 (0.7)   | -1.21 (0.92) | 2.26 (1.25) |
| Diaconis-Ylvisaker                      | -7.58 (5.66) | -0.06 (0.05) | 0.62 (0.55) | -0.51 (0.56) | -1.47 (0.68) | -1.09 (0.9)  | 2.03 (1.22) |
| Clogg et al. (1991)                     | -7.67 (5.7)  | -0.06 (0.05) | 0.63 (0.55) | -0.52 (0.57) | -1.47 (0.68) | -1.1 (0.9)   | 2.06 (1.22) |
| Firth (1993)                            | -7.4 (5.66)  | -0.06 (0.05) | 0.62 (0.55) | -0.53 (0.56) | -1.45 (0.68) | -1.1 (0.9)   | 2 (1.22)    |
| Kenne Pagui, Salvan, and Sartori (2017) | -7.64 (5.72) | -0.06 (0.05) | 0.64 (0.56) | -0.54 (0.57) | -1.48 (0.68) | -1.13 (0.91) | 2.06 (1.23) |

MLE and BR estimates

We assess the properties of different approaches relying on a simulation
study. Specifically, 10000 artificial datasets are generated from the
maximum-likelihood estimates. The approaches illustrated on Table 1 are
estimated for each replications, and we evaluate their performance in
terms of bias and root mean squared error. Note that the simulation
requires roughly 3 minutes on a 2020 Macbook Pro with M1 processor
(`aarch64-apple-darwin20`) running R 4.1.1 linked with `openblas`.

For convenience, we store the results in the file
[`birthweight_sim.RData`](https://github.com/tommasorigon/logistic-bias-reduction/blob/main/BIRTHWEIGHT/birthweight_sim.RData),
which can be loaded directly; refer to
[birthweight.Rmd](https://raw.githubusercontent.com/tommasorigon/logistic-bias-reduction/main/BIRTHWEIGHT/birthweight.Rmd),
for a complete illustration of the code adapted from Kosmidis, Kenne
Pagui, and Sartori (2020).

``` r
load("birthweight_sim.RData")
```

## Bias

|                                         |       |       |       |       |       |       |      |
|:----------------------------------------|------:|------:|------:|------:|------:|------:|-----:|
| MLE                                     | -1.42 | -0.01 |  0.09 | -0.04 | -0.18 | -0.12 | 0.34 |
| Diaconis-Ylvisaker                      | -0.08 |  0.00 | -0.01 |  0.03 | -0.01 |  0.03 | 0.00 |
| Clogg et al. (1991)                     | -0.22 |  0.00 |  0.00 |  0.02 | -0.01 |  0.03 | 0.05 |
| Firth (1993)                            | -0.08 |  0.00 |  0.01 |  0.00 |  0.00 |  0.00 | 0.02 |
| Kenne Pagui, Salvan, and Sartori (2017) | -0.38 |  0.00 |  0.03 | -0.01 | -0.06 | -0.03 | 0.10 |

Bias

## RMSE

|                                         |      |      |      |      |      |      |      |
|:----------------------------------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| MLE                                     | 6.88 | 0.06 | 0.64 | 0.65 | 0.81 | 1.12 | 1.49 |
| Diaconis-Ylvisaker                      | 5.71 | 0.05 | 0.54 | 0.56 | 0.71 | 0.96 | 1.22 |
| Clogg et al. (1991)                     | 5.83 | 0.05 | 0.55 | 0.57 | 0.70 | 0.97 | 1.25 |
| Firth (1993)                            | 5.94 | 0.05 | 0.57 | 0.58 | 0.71 | 0.95 | 1.28 |
| Kenne Pagui, Salvan, and Sartori (2017) | 6.12 | 0.06 | 0.58 | 0.60 | 0.78 | 1.01 | 1.31 |

RMSE

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Clogg1991" class="csl-entry">

Clogg, Clifford C., Donald B. Rubin, Nathaniel Schenker, Bradley
Schultz, and Lynn Weidman. 1991. “<span class="nocase">Multiple
imputation of industry and occupation codes in census public-use samples
using Bayesian logistic regression</span>.” *J. Am. Statist. Assoc.* 86
(413): 68–78. <https://doi.org/10.1080/01621459.1991.10475005>.

</div>

<div id="ref-Firth1993" class="csl-entry">

Firth, David. 1993. “<span class="nocase">Bias reduction of maximum
likelihood estimates</span>.” *Biometrika* 80 (1): 27–38.

</div>

<div id="ref-Pagui2017" class="csl-entry">

Kenne Pagui, E. C., A. Salvan, and N. Sartori. 2017. “<span
class="nocase">Median bias reduction of maximum likelihood
estimates</span>.” *Biometrika* 104 (4): 923–38.

</div>

<div id="ref-Kosmidis2021" class="csl-entry">

Kosmidis, Ioannis, and David Firth. 2021. “<span
class="nocase">Jeffreys-prior penalty, finiteness and shrinkage in
binomial-response generalized linear models</span>.” *Biometrika* 108
(1): 71–82. <https://doi.org/10.1093/biomet/asaa052>.

</div>

<div id="ref-Kosmidis2020" class="csl-entry">

Kosmidis, Ioannis, Euloge Clovis Kenne Pagui, and Nicola Sartori. 2020.
“<span class="nocase">Mean and median bias reduction in generalized
linear models</span>.” *Statist. Comp.* 30 (1): 43–59.
<https://doi.org/10.1007/s11222-019-09860-6>.

</div>

</div>
