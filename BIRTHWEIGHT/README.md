# Birthweight Study

This tutorial is devoted to replicating simulations from Kosmidis, Kenne
Pagui, and Sartori (2020), available in Table 2 of the article Data
comprises *n* = 100 births and the binary outcome of interest is a
dichotomization of infant birthweight (low against normal). A logistic
regression is employed with covariates pertaining the mother, namely:
age, race, smoking habitudes, history of premature labour and
hypertension and the logarithm of the weight at birth of the mother.

Similarly to the `endometrial` [case study](../ENDOMETRIAL/), we compare
the DY regularizaion with Firth (1993), Kenne Pagui, Salvan, and Sartori
(2017) and Clogg et al. (1991), relying on the package `brglm2`
(Kosmidis and Firth 2021).

Estimated regression coefficients are reported in the following table.
Note that maximum-likelihood estimate exists finite in this example.

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

|                                         | (Intercept) |    age | racewhite |  smoke |    ptl |     ht | loglwt |
|:----------------------------------------|------------:|-------:|----------:|-------:|-------:|-------:|-------:|
| MLE                                     |      -8.496 | -0.067 |     0.690 | -0.560 | -1.603 | -1.211 |  2.262 |
| Diaconis-Ylvisaker                      |      -7.584 | -0.060 |     0.618 | -0.505 | -1.469 | -1.091 |  2.026 |
| Firth (1993)                            |      -7.401 | -0.061 |     0.622 | -0.531 | -1.446 | -1.104 |  1.998 |
| Kenne Pagui, Salvan, and Sartori (2017) |      -7.641 | -0.062 |     0.638 | -0.538 | -1.481 | -1.134 |  2.059 |
| Clogg et al. (1991)                     |      -7.665 | -0.061 |     0.629 | -0.515 | -1.466 | -1.098 |  2.057 |

MLE and BR estimates

We assess properties of estimates relying on a simulation study. 10000
artificial dataset are generated from the maximum-likelihood estimates.
Note that the simulation require rougly 3 minutes on a 2020 Macbook Pro
with M1 processor (`aarch64-apple-darwin20`) running R 4.1.1 linked with
`openblas`. For convenience, we store the results in the file
`birthweight_sim.RData`, provided within this folder.

## Bias

|                                         |         |       |       |       |        |        |       |
|:----------------------------------------|--------:|------:|------:|------:|-------:|-------:|------:|
| MLE                                     | -141.56 | -0.79 |  8.90 | -3.61 | -18.04 | -12.20 | 34.13 |
| Diaconis-Ylvisaker                      |   -8.28 | -0.13 |  0.75 | -0.39 |  -0.02 |  -0.04 |  2.19 |
| Firth (1993)                            |  -38.16 | -0.21 |  2.79 | -1.20 |  -5.57 |  -3.27 |  9.63 |
| Kenne Pagui, Salvan, and Sartori (2017) |  -21.67 | -0.04 |  0.39 |  1.93 |  -0.78 |   2.92 |  4.85 |
| Clogg et al. (1991)                     |   -8.13 |  0.07 | -1.09 |  3.04 |  -1.01 |   3.45 |  0.44 |

Bias x 100

## Rmse

|                                         |      |      |      |      |      |      |      |
|:----------------------------------------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|
| MLE                                     | 6.88 | 0.06 | 0.64 | 0.65 | 0.81 | 1.12 | 1.49 |
| Diaconis-Ylvisaker                      | 5.94 | 0.05 | 0.57 | 0.58 | 0.71 | 0.95 | 1.28 |
| Firth (1993)                            | 6.12 | 0.06 | 0.58 | 0.60 | 0.78 | 1.01 | 1.31 |
| Kenne Pagui, Salvan, and Sartori (2017) | 5.83 | 0.05 | 0.55 | 0.57 | 0.70 | 0.97 | 1.25 |
| Clogg et al. (1991)                     | 5.71 | 0.05 | 0.54 | 0.56 | 0.71 | 0.96 | 1.22 |

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
