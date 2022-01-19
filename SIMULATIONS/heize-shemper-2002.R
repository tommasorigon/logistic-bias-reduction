## Reproduce Section 3 from heinz - Schemper (2002, stat in medicine)
rm(list = ls())
N = c(30,50,100)
P = c(3,5,10)
OR = c(1,2,4,16)
pr_X = c(0.5,0.25)
pr_Y = c(0.5, 0.25)
alpha_grid = seq(-4,4,l=100)

out = list()
for(id_n in 1:3){
  for(id_p in 1:3){
    for(id_or in 1:3){
      for(id_X in 1:2){
        for(id_Y in 1:2){


          set.seed(1991)
          n = N[id_n]
          p = P[id_p]
          or = OR[id_or]
          b = rep(log(or),p)
          X = matrix(rbinom(n*p,1,pr_X[id_X]),n)
          y_bar = sapply(alpha_grid, function(x) mean(plogis(X %*% b + x)))
          b0 = alpha_grid[which.min(abs(y_bar - pr_Y[id_Y]))]

          X = cbind(1,X)
          beta = c(b0,b)
          Nsim = 5e3
          ml <- br <- clogg <- dy <- matrix(NA, nrow = Nsim, ncol = length(beta))

          # Replicate 1000 times
          for(i in 1:Nsim){
            if (i %% 10 == 0) print(i)
            y = rbinom(n,1,prob = plogis(X %*% beta))
            mle <- glm(y ~ X - 1, family = binomial("logit"))
            m = NROW(model.matrix(mle))
            p = NCOL(model.matrix(mle))

            y_dy <- p / (p + m) * 0.5 + m / (p + m) * y
            y_clogg <- p / (p + m) * mean(y) + m / (p + m) * y

            # DY ESTIMATE
            fit_fast_dy <- fastLR(X, y_dy)
            fit_fast_clogg <- fastLR(X, y_clogg)

            # FIRTH (1993)
            fit_firth <- glm(y ~ -1 + X,
                             family = binomial("logit"),
                             method = "brglmFit", type = "AS_mean"
            )



            ml[i, ] <- coef(mle)
            dy[i, ] <- fit_fast_dy$coefficients
            clogg[i, ] <- fit_fast_clogg$coefficients
            br[i, ] <- coef(fit_firth)
          }


          ## BIAS in beta parameterization
          bias.beta <- data.frame(
            ml = colMeans(ml[]),
            br = colMeans(br),
            clogg = colMeans(clogg),
            dy = colMeans(dy)
          ) - beta

          ## MSE in beta parameterization
          mse.beta <- data.frame(
            ml = rowMeans((t(ml[, ]) - beta)^2),
            br = rowMeans((t(br) - beta)^2),
            clogg = rowMeans((t(clogg) - beta)^2),
            dy = rowMeans((t(dy) - beta)^2)
          )

          # save results
          settings = list("n" = n,"p" = p, "or" = or, "Bx" = pr_X[id_X], "By" = pr_Y[id_Y])
          out[[length(out) + 1]] = list("settings" = settings, "bias" = bias.beta, "mse" = mse.beta)
        }
      }
    }
  }
}






df = reshape2::melt(list("Bias" = bias.beta, "RMSE" = sqrt(mse.beta)))

df$variable = factor(df$variable, labels = c("MLE", "Firth (1993)", "Clogg (1991)", "DY"))
library(ggplot2)
pl = ggplot(df) +
  geom_boxplot(aes(y = value, x = variable,  group = variable), color = "gray50") +
  geom_jitter(aes(y = value, x = variable,  group = variable), alpha = .2, color = "gray") +
  facet_wrap(~L1, scales = "free") +
  theme_bw() + xlab("") + ylab("") + theme(legend.position = "none")
pl
