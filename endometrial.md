## Endometrial cancer study

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

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

    beta_clogg <- coef(fit_clogg)
    sd_clogg  <- summary(fit_clogg)$coefficients[,2]


    # DIACONIS & YLVISAKER
    fit_dy <- glm(HG_DY ~ NV + PI + EH, data = endometrial, family = binomial("logit"))

    ## Warning in eval(family$initialize): non-integer #successes in a binomial glm!

    beta_dy <- coef(fit_dy)
    sd_dy  <- summary(fit_dy)$coefficients[,2]

    estimates <- rbind(beta_mle, beta_firth, beta_kp, beta_clogg, beta_dy)
    std_errors <- rbind(sd_mle, sd_firth, sd_kp, sd_clogg, sd_dy)
    kable(estimates, digits = 3)

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">(Intercept)</th>
<th style="text-align: right;">NV</th>
<th style="text-align: right;">PI</th>
<th style="text-align: right;">EH</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">beta_mle</td>
<td style="text-align: right;">4.305</td>
<td style="text-align: right;">Inf</td>
<td style="text-align: right;">-0.042</td>
<td style="text-align: right;">-2.903</td>
</tr>
<tr class="even">
<td style="text-align: left;">beta_firth</td>
<td style="text-align: right;">3.775</td>
<td style="text-align: right;">2.929</td>
<td style="text-align: right;">-0.035</td>
<td style="text-align: right;">-2.604</td>
</tr>
<tr class="odd">
<td style="text-align: left;">beta_kp</td>
<td style="text-align: right;">3.969</td>
<td style="text-align: right;">3.869</td>
<td style="text-align: right;">-0.039</td>
<td style="text-align: right;">-2.708</td>
</tr>
<tr class="even">
<td style="text-align: left;">beta_clogg</td>
<td style="text-align: right;">3.622</td>
<td style="text-align: right;">3.223</td>
<td style="text-align: right;">-0.034</td>
<td style="text-align: right;">-2.511</td>
</tr>
<tr class="odd">
<td style="text-align: left;">beta_dy</td>
<td style="text-align: right;">3.579</td>
<td style="text-align: right;">3.431</td>
<td style="text-align: right;">-0.034</td>
<td style="text-align: right;">-2.458</td>
</tr>
</tbody>
</table>

    kable(std_errors, digits = 3)

<table>
<thead>
<tr class="header">
<th style="text-align: left;"></th>
<th style="text-align: right;">(Intercept)</th>
<th style="text-align: right;">NV</th>
<th style="text-align: right;">PI</th>
<th style="text-align: right;">EH</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">sd_mle</td>
<td style="text-align: right;">1.637</td>
<td style="text-align: right;">Inf</td>
<td style="text-align: right;">0.044</td>
<td style="text-align: right;">0.846</td>
</tr>
<tr class="even">
<td style="text-align: left;">sd_firth</td>
<td style="text-align: right;">1.489</td>
<td style="text-align: right;">1.551</td>
<td style="text-align: right;">0.040</td>
<td style="text-align: right;">0.776</td>
</tr>
<tr class="odd">
<td style="text-align: left;">sd_kp</td>
<td style="text-align: right;">1.552</td>
<td style="text-align: right;">2.298</td>
<td style="text-align: right;">0.042</td>
<td style="text-align: right;">0.803</td>
</tr>
<tr class="even">
<td style="text-align: left;">sd_clogg</td>
<td style="text-align: right;">1.471</td>
<td style="text-align: right;">1.722</td>
<td style="text-align: right;">0.040</td>
<td style="text-align: right;">0.761</td>
</tr>
<tr class="odd">
<td style="text-align: left;">sd_dy</td>
<td style="text-align: right;">1.459</td>
<td style="text-align: right;">1.893</td>
<td style="text-align: right;">0.040</td>
<td style="text-align: right;">0.748</td>
</tr>
</tbody>
</table>
