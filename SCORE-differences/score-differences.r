---
title: "Endometrial cancer study"
author: "Tommaso Rigon and Emanuele Aliverti"
bibliography: ["../biblio.bib"]
output: html_document
---

rm(list = ls())
library(knitr)
library(jtools)
library(brglm2)
library(enrichwith)
data(endometrial)


y <- endometrial$HG
X <- model.matrix(HG ~ NV + PI + EH, data = endometrial)
X[,3] = (X[,3] - min(X[,3]))/  (max(X[,3]) - min(X[,3]))
X[,4] = (X[,4] - min(X[,4]))/  (max(X[,4]) - min(X[,4]))
p <- ncol(X)
m <- nrow(X)
endometrial$HG_DY <- p / (p + m) * 0.5 + m / (p + m) * endometrial$HG
endometrial$PIs = X[,3]
endometrial$EHs = X[,4]

fit_glm =  glm(HG ~ NV + PIs + EHs, family = binomial("logit"), data = endometrial)
fit_dy <- glm(HG_DY ~ NV + PIs + EHs, family = binomial("logit"), data = endometrial)
fit_firth <- glm(HG ~NV + PIs + EHs, family = binomial("logit"), method = "brglmFit", type = "AS_mean", data = endometrial)

score_glm = function(beta,X = model.matrix(fit_glm), y = fit_glm$y){
   pr = c(plogis(tcrossprod(beta,X)))
   colSums((y-pr)*X)
}

score_dy = function(beta,X = model.matrix(fit_glm), y = fit_glm$y) {
 # build probabilities 
  pr = c(plogis(tcrossprod(beta,X)))
  score_glm(beta,X,y) - colSums( NCOL(X)/NROW(X) * (pr-0.5) * X) 
}

# crea correzione firth a partire dal glm
score_firth = function(beta, X = model.matrix(fit_glm),y =fit_glm$y) {
 # build probabilities 
  pr = c(plogis(tcrossprod(beta,X)))
  ppr = pr * (1-pr)
 # build leverages 
  M = solve(crossprod(X*ppr, X))
  S = X*sqrt(ppr)
  W = S %*% M %*% t(S)
  h = diag(W)
  
  bias = h*(pr-0.5)*X
  score_glm(beta, X,y) - colSums(bias)
  
}

score_glm(coef(fit_glm))
score_firth(coef(fit_dy))
score_dy(coef(fit_firth))
c2 = coef(fit_dy)
c2[2] = 3
score_dy(c2)
## MLE (add lines)
estimates = lapply(list(DY=fit_dy,FI=fit_firth,ML=fit_glm), function(x) data.frame(t(coef(x))))

# get a plausible range of values
par_ranges = apply(confint(fit_firth, level = 0.5), 1, function(x) seq(x[1], x[2], l = 100))
par
# fix the other components at their maximum
plot(hatvalues(fit_dy), hatvalues(fit_firth))
plot(hatvalues(fit_dy), hatvalues(fit_glm))
plot(fitted(fit_dy), fitted(fit_firth))

abline(a = 0,b=1,col=2)

hatvalues(fit_dy)

# other components should be fixed at their maximum
sc_fi

sc_fi = t(apply(par_ranges, 1, function(x) score_firth(x)))
sc_dy = t(apply(par_ranges, 1, function(x) score_dy(x)))
sc_ml = t(apply(par_ranges, 1, function(x) score_glm(x)))


apply(abs(sc_fi - sc_dy),2, max)

tmp = reshape2::melt(par_ranges)
tmp$FI = reshape2::melt(sc_fi)$value
tmp$DY = reshape2::melt(sc_dy)$value
tmp$ml = reshape2::melt(sc_ml)$value
tmp
## MLE (add lines)
estimates = lapply(list(DY=fit_dy,FI=fit_firth,ML=fit_glm), function(x) data.frame(t(coef(x))))
est_plot = reshape2::melt(estimates)
est_plot$variable = factor(est_plot$variable, labels = levels(tmp$Var2))
est_plot$Var2 = est_plot$variable
# rimuovi MLE infinito
est_plot$value[est_plot$value > 10] = NA

require(ggplot2)
p0 = ggplot(tmp) + 
  geom_line(aes(x = value, y = FI), color = "green") +
  geom_line(aes(x = value, y = DY), color = "red") +
  geom_line(aes(x = value, y = ml), color = "blue") +
 # geom_vline(data = est_plot, aes(xintercept = value, color = L1), lty = 2)+
  facet_wrap(~Var2, scales = "free") +
  geom_hline(yintercept = 0, lty = 3)+
  theme_minimal()
  # plot the score functions
p0
ggsave(p0, file = "/tmp/pl2.png")



