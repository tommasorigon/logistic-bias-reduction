rm(list = ls())

## Compute standard errors
load("full_sim.RData")
load("sur-candes.RData")
se_compute = function(beta,X){
  pi = plogis(X %*% beta)
  W = pi * (1-pi)
  se = sqrt(diag(solve(crossprod(c(W) * X,X))))
  return(se)
}

ml_se = apply(ml, 1, function(x) se_compute(x,X))
dy_se = apply(dy, 1, function(x) se_compute(x,X))
br_se = apply(br, 1, function(x) se_compute(x,X))
clogg_se = apply(clogg, 1, function(x) se_compute(x,X))


## Coverage of 95% Wald confidence intervals
ml.ci.l <- ml - qnorm(0.975) * t(ml_se)
ml.ci.u <- ml + qnorm(0.975) * t(ml_se)

br.ci.l <- br - qnorm(0.975) * t(br_se)
br.ci.u <- br + qnorm(0.975) * t(br_se)

dy.ci.l <- dy - qnorm(0.975) * t(dy_se)
dy.ci.u <- dy + qnorm(0.975) * t(dy_se)

clogg.ci.l <- clogg - qnorm(0.975) * t(clogg_se)
clogg.ci.u <- clogg + qnorm(0.975) * t(clogg_se)


coverage <- data.frame(ml = rowMeans(t(ml.ci.l) < beta &  t(ml.ci.u) > beta),
                       br = rowMeans(t(br.ci.l) < beta &  t(br.ci.u) > beta),
                       dy = rowMeans(t(dy.ci.l) < beta &  t(dy.ci.u) > beta),
                       clogg = rowMeans(t(clogg.ci.l) < beta &  t(clogg.ci.u) > beta))

## probability of underestimation
PU <- data.frame(ml = rowMeans(t(ml) < beta),
                 br = rowMeans(t(br) < beta),
                 dy = rowMeans(t(dy) < beta),
                 clogg = rowMeans(t(clogg) < beta))


coverage_p = round(coverage*100,1)
pu_p = round(PU*100,1)
rmse.beta = sqrt(mse.beta)
pu_p$x = coverage_p$x = rmse.beta$x = bias.beta$x = mae.beta$x = 1:NROW(bias.beta)

df = reshape2::melt(list("Bias" = bias.beta,
                         "RMSE" = rmse.beta,
                         "MAE" = mae.beta,
                         "wald" = coverage_p), id.var = "x")
df$variable = factor(df$variable, levels = c("ml", "dy", "clogg", "br"),
                     labels=                     c("MLE","DY",  "Clogg (1991)", "Firth (1993)"),ord=T)
(df$variable)
df$variable = factor(df$variable, labels = c("Maximum Likelihood","Diaconis-Ylvisaker",  "Clogg (1991)", "Firth (1993)") )
df$b = cut(df$x,breaks = c(0,25,50,200),labels = c("10", "-10", "0"))
library(dplyr)
ff = function(x) {
  tmp = round(quantile(x, c(.25,.5,.75)),2)
  # sprintf("%g (%g,%g)", tmp[2],tmp[1],tmp[3])
  sprintf("%.2f (%.2f)", tmp[2],abs(tmp[3] -tmp[1]))

}
res = df %>% group_by(L1,b,variable) %>%
  summarize(x=ff(value))
tab_cleaned = reshape2::dcast(b + variable ~ L1,data = res)
tab_cleaned
xtable::print.xtable(xtable::xtable(tab_cleaned),row.names = F)



