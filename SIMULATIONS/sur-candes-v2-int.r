rm(list = ls())
set.seed(11235)
n <- 1000
p <- 200
b0 = c(1, -1, 2,-2, 0)
beta = rep(b0,each = p/5)
X <- matrix(rnorm(n * (p), 0, sqrt(1 / n)), n, p)
X = cbind(1,X[,-1])
y <- rbinom(n, 1, plogis(X %*% beta))
mean(y)
m0 <- glm(y ~ X - 1, family = binomial("logit"))
Nsim <- 5e3
y_sim <- matrix(NA,n,Nsim)
# generate from the true values
y_sim <- apply(y_sim,2, function(x) rbinom(n, 1, plogis(X %*% beta)))


ml <- br <- clogg <- dy <- matrix(NA, nrow = Nsim, ncol = length(beta))

library(brglm2)
library(RcppNumerical)
m <- nrow(X)
for (i in 1:Nsim) {
  if (i %% 5 == 0) print(i)
  y <- y_sim[, i]

  mle <- fastLR(X,y)

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



  ml[i, ] <- mle$coefficients
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

## MAE in beta parameterization
rowMedian = function(x) apply(x,1,median)
mae.beta <- data.frame(
  ml = rowMedian(abs(t(ml[, ]) - beta)),
  br = rowMedian(abs(t(br) - beta)),
  clogg = rowMedian(abs(t(clogg) - beta)),
  dy = rowMedian(abs(t(dy) - beta))
)

## OR
## BIAS in OR parameterization
bias.or <- data.frame(
  ml = colMeans(exp(ml[])),
  br = colMeans(exp(br)),
  clogg = colMeans(exp(clogg)),
  dy = colMeans(exp(dy))
) - exp(beta)


or = exp(beta)

## MSE in OR parameterization
mse.or <- data.frame(
  ml = rowMeans((t(exp(ml[, ])) -or)^2),
  br = rowMeans((t(exp(br)) -or)^2),
  clogg = rowMeans((t(exp(clogg)) -or)^2),
  dy = rowMeans((t(exp(dy)) -or)^2)
)



## MAE in OR parameterization
mae.or <- data.frame(
  ml = rowMedian(abs(t(exp(ml[, ])) -or)),
  br = rowMedian(abs(t(exp(br)) -or)),
  clogg = rowMedian(abs(t(exp(clogg)) -or)),
  dy = rowMedian(abs(t(exp(dy)) -or))
)



beta_res = c("mse.beta", "bias.beta", "mae.beta")
or_res = c("mse.or", "bias.or", "mae.or")


save(list = c(beta_res, or_res), file = "sur-candes-v2-int.RData")
save.image(file = "full_sim-v2-int.RData")


