## Endometrial cancer study

``` r
library(brglm2)
library(knitr)

rm(list = ls())
data(endometrial)

# DATA
y <- endometrial$HG
X <- model.matrix(HG ~ NV + PI + EH, data = endometrial)
p <- ncol(X)
m <- nrow(X)
endometrial$HG_CLOGG <- p / (p + m) * mean(endometrial$HG) + m / (p + m) * endometrial$HG
endometrial$HG_DY <- p / (p + m) * 0.5 + m / (p + m) * endometrial$HG

# MLE ESTIMATE
fit_mle <- glm(HG ~ NV + PI + EH, data = endometrial, family = binomial("logit"))
beta_mle <- coef(fit_mle)
sd_mle  <- summary(fit_mle)$coefficients[,2]

# I manually set these parameters to + Inf 
beta_mle[2] <- sd_mle[2] <- Inf

# FIRTH (1993)
fit_firth <- glm(HG ~ NV + PI + EH,
  data = endometrial, family = binomial("logit"),
  method = "brglmFit", type = "AS_mean"
)
beta_firth <- coef(fit_firth)
sd_firth  <- summary(fit_firth)$coefficients[,2]

# KENNE PAGUI ET AL. (2017)
fit_kp <- glm(HG ~ NV + PI + EH,
  data = endometrial, family = binomial("logit"),
  method = "brglmFit", type = "AS_median"
)
beta_kp <- coef(fit_kp)
sd_kp  <- summary(fit_kp)$coefficients[,2]

# CLOGG ET AL. (1991)
fit_clogg <- glm(HG_CLOGG ~ NV + PI + EH, data = endometrial, family = binomial("logit"))
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
beta_clogg <- coef(fit_clogg)
sd_clogg  <- summary(fit_clogg)$coefficients[,2]


# DIACONIS & YLVISAKER
fit_dy <- glm(HG_DY ~ NV + PI + EH, data = endometrial, family = binomial("logit"))
```

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

``` r
beta_dy <- coef(fit_dy)
sd_dy  <- summary(fit_dy)$coefficients[,2]

estimates <- rbind(beta_mle, beta_firth, beta_kp, beta_clogg, beta_dy)
std_errors <- rbind(sd_mle, sd_firth, sd_kp, sd_clogg, sd_dy)
kable(estimates, digits = 3)
```

|            | (Intercept) |    NV |     PI |     EH |
|:-----------|------------:|------:|-------:|-------:|
| beta_mle   |       4.305 |   Inf | -0.042 | -2.903 |
| beta_firth |       3.775 | 2.929 | -0.035 | -2.604 |
| beta_kp    |       3.969 | 3.869 | -0.039 | -2.708 |
| beta_clogg |       3.622 | 3.223 | -0.034 | -2.511 |
| beta_dy    |       3.579 | 3.431 | -0.034 | -2.458 |

``` r
kable(std_errors, digits = 3)
```

|          | (Intercept) |    NV |    PI |    EH |
|:---------|------------:|------:|------:|------:|
| sd_mle   |       1.637 |   Inf | 0.044 | 0.846 |
| sd_firth |       1.489 | 1.551 | 0.040 | 0.776 |
| sd_kp    |       1.552 | 2.298 | 0.042 | 0.803 |
| sd_clogg |       1.471 | 1.722 | 0.040 | 0.761 |
| sd_dy    |       1.459 | 1.893 | 0.040 | 0.748 |
