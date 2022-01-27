## Endometrial cancer study

*(This tutorial illustrates main results of the analysis in an
easy-to-read format. Complete code associated with this page is
available in the file*
[endometrial.Rmd](https://github.com/tommasorigon/logistic-bias-reduction/blob/main/ENDOMETRIAL/endometrial.Rmd)),

The endometrial cancer study is illustrated in Heinze and Schemper
(2002), and focuses on *n* = 79 patients to evaluate the relationship
between the histology of the endometrium (low against high), and three
risk factors: neovasculation, pulsatility index of arteria uterina, and
endometrium height. As shown in Heinze and Schemper (2002), maximum
likelihood for the coefficient *β*<sub>2</sub> (neovasculation) does not
exist.

``` r
library(brglm2)
data(endometrial)
fit_mle <- glm(HG ~ NV + PI + EH, data = endometrial, family = binomial("logit"))
```

    ## Standard errors: MLE
    ## ---------------------------------------------------
    ##                      Est.      S.E.   z val.      p
    ## ----------------- ------- --------- -------- ------
    ## (Intercept)          4.30      1.64     2.63   0.01
    ## NV                  18.19   1715.75     0.01   0.99
    ## PI                  -0.04      0.04    -0.95   0.34
    ## EH                  -2.90      0.85    -3.43   0.00
    ## ---------------------------------------------------

## Regularization and bias correction

As illustrated in Section 1 of the article, our approach relies on a
Diaconis-Ylvisaker conjugate prior with carefully chosen
hyperparameters. The resulting penalized likelihood can be expressed as
a genuine Binomial likelihood, obtained replacing the Binomial
observations via pseudo-counts. We construct such quantities in the
variable `HG_DY`.

``` r
y <- endometrial$HG
X <- model.matrix(HG ~ NV + PI + EH, data = endometrial)
p <- ncol(X)
m <- nrow(X)
endometrial$HG_DY <- p / (p + m) * 0.5 + m / (p + m) * endometrial$HG
```

Our method can be estimated directly via `glm`, or other more efficient
implementations for logistic regression models (refer to
[HIGH-DIMENSIONAL-SYNTHETIC](../HIGH-DIMENSIONAL-SYNTHETIC) for further
examples in higher dimensions). Note that pseudo counts are not
necessary integers, and this might prompt a warning when `glm` is called
with `binomial` family. Such warnings can be avoided relying on the
`quasi` family with link and variance functions corresponding to
binomial likelihood. In both cases, penalized likelihood optimization
provide finite values for all coefficients, in virtue of Theorem 1 of
the article.

``` r
fit_dy <- glm(HG_DY ~ NV + PI + EH, data = endometrial, family = binomial("logit"))
# Avoid warnings
# fit_dy_now <- glm(HG_DY ~ NV + PI + EH, data = endometrial, family = quasi(link = "logit",variance = "mu(1-mu)"))
```

    ## Standard errors: MLE
    ## ------------------------------------------------
    ##                      Est.   S.E.   z val.      p
    ## ----------------- ------- ------ -------- ------
    ## (Intercept)          3.58   1.46     2.45   0.01
    ## NV                   3.43   1.89     1.81   0.07
    ## PI                  -0.03   0.04    -0.86   0.39
    ## EH                  -2.46   0.75    -3.29   0.00
    ## ------------------------------------------------

Table 1 in the paper compares the proposed approach with state-of-the
art methods for bias-reduction, such as Firth (1993) and Kenne Pagui,
Salvan, and Sartori (2017), conveniently implemented in the R package
`brglm2` (Kosmidis 2021); refer also to Kosmidis, Kenne Pagui, and
Sartori (2020). We also compare the proposed method with the approach of
(Clogg et al. 1991), which relies on a conjugate prior specification
inducing pseudo-counts and shrinking estimates toward the empirical
proportion.

``` r
# CLOGG ET AL. (1991)
endometrial$HG_CLOGG <- p / (p + m) * mean(endometrial$HG) + m / (p + m) * endometrial$HG
fit_clogg <- glm(HG_CLOGG ~ NV + PI + EH, data = endometrial, family = binomial("logit"))

# FIRTH (1993)
fit_firth <- glm(HG ~ NV + PI + EH,
  data = endometrial, family = binomial("logit"),
  method = "brglmFit", type = "AS_mean"
)

# KENNE PAGUI ET AL. (2017)
fit_kp <- glm(HG ~ NV + PI + EH,
  data = endometrial, family = binomial("logit"),
  method = "brglmFit", type = "AS_median"
)
```

|                                         |               |               |                |                |
|:----------------------------------------|:--------------|:--------------|:---------------|:---------------|
| MLE                                     | 4.305 (1.637) | Inf (Inf)     | -0.042 (0.044) | -2.903 (0.846) |
| Diaconis-Ylvisaker                      | 3.579 (1.459) | 3.431 (1.893) | -0.034 (0.04)  | -2.458 (0.748) |
| Clogg et al. (1991)                     | 3.622 (1.471) | 3.223 (1.722) | -0.034 (0.04)  | -2.511 (0.761) |
| Firth (1993)                            | 3.775 (1.489) | 2.929 (1.551) | -0.035 (0.04)  | -2.604 (0.776) |
| Kenne Pagui, Salvan, and Sartori (2017) | 3.969 (1.552) | 3.869 (2.298) | -0.039 (0.042) | -2.708 (0.803) |

Table 1 of the paper: Estimated regression coefficients on the
endometrial cancer study (and standard errors)

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

<div id="ref-Heinze2002" class="csl-entry">

Heinze, Georg, and Michael Schemper. 2002. “<span class="nocase">A
solution to the problem of separation in logistic regression</span>.”
*Statist. Med.* 21 (16): 2409–19. <https://doi.org/10.1002/sim.1047>.

</div>

<div id="ref-Pagui2017" class="csl-entry">

Kenne Pagui, E. C., A. Salvan, and N. Sartori. 2017. “<span
class="nocase">Median bias reduction of maximum likelihood
estimates</span>.” *Biometrika* 104 (4): 923–38.

</div>

<div id="ref-brglm2" class="csl-entry">

Kosmidis, Ioannis. 2021. *<span class="nocase">brglm2</span>: Bias
Reduction in Generalized Linear Models*.
<https://CRAN.R-project.org/package=brglm2>.

</div>

<div id="ref-Kosmidis2020" class="csl-entry">

Kosmidis, Ioannis, Euloge Clovis Kenne Pagui, and Nicola Sartori. 2020.
“<span class="nocase">Mean and median bias reduction in generalized
linear models</span>.” *Statist. Comp.* 30 (1): 43–59.
<https://doi.org/10.1007/s11222-019-09860-6>.

</div>

</div>
