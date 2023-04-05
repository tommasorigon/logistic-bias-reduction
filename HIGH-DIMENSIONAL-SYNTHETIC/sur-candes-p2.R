rm(list = ls())
# This version focuses on a ratio p/n closer to 1
set.seed(1991)
n <- 1000
p <- 300
b0 = c(3,1.5,0, -1.5,-3)
beta = rep(b0,each = p/5)
X <- matrix(rnorm(n * (p), 0, sqrt(1 / n)), n, p)
y <- rbinom(n, 1, plogis(X %*% beta))
m0 <- glm(y ~ X - 1, family = binomial("logit"))
Nsim <- 5e3
y_sim <- matrix(NA,n,Nsim)
# generate from the true values
y_sim <- apply(y_sim,2, function(x) rbinom(n, 1, plogis(X %*% beta)))


ml <- br <- dy <- matrix(NA, nrow = Nsim, ncol = length(beta))

library(brglm2)
library(RcppNumerical)
m <- nrow(X)
for (i in 1:Nsim) {
  if (i %% 5 == 0) print(i)
  y <- y_sim[, i]

  mle <- fastLR(X,y)

  y_dy <- p / (p + m) * 0.5 + m / (p + m) * y

  # DY ESTIMATE
  fit_fast_dy <- fastLR(X, y_dy)

  # FIRTH (1993)
  fit_firth <- glm(y ~ -1 + X,
    family = binomial("logit"),
    method = "brglmFit", type = "AS_mean",
    control = list(maxit = 100)
  )



  ml[i, ] <- mle$coefficients
  dy[i, ] <- fit_fast_dy$coefficients
  br[i, ] <- coef(fit_firth)
}


## BIAS in beta parameterization
bias.beta <- data.frame(
  ml = colMeans(ml[]),
  br = colMeans(br),
  dy = colMeans(dy)
) - beta

## MSE in beta parameterization
mse.beta <- data.frame(
  ml = rowMeans((t(ml[, ]) - beta)^2),
  br = rowMeans((t(br) - beta)^2),
  dy = rowMeans((t(dy) - beta)^2)
)

## MAE in beta parameterization
rowMedian = function(x) apply(x,1,median)
mae.beta <- data.frame(
  ml = rowMedian(abs(t(ml[, ]) - beta)),
  br = rowMedian(abs(t(br) - beta)),
  dy = rowMedian(abs(t(dy) - beta))
)

## OR
## BIAS in OR parameterization
bias.or <- data.frame(
  ml = colMeans(exp(ml[])),
  br = colMeans(exp(br)),
  dy = colMeans(exp(dy))
) - exp(beta)


or = exp(beta)

## MSE in OR parameterization
mse.or <- data.frame(
  ml = rowMeans((t(exp(ml[, ])) -or)^2),
  br = rowMeans((t(exp(br)) -or)^2),
  dy = rowMeans((t(exp(dy)) -or)^2)
)



## MAE in OR parameterization
mae.or <- data.frame(
  ml = rowMedian(abs(t(exp(ml[, ])) -or)),
  br = rowMedian(abs(t(exp(br)) -or)),
  dy = rowMedian(abs(t(exp(dy)) -or))
)



beta_res = c("mse.beta", "bias.beta", "mae.beta")
or_res = c("mse.or", "bias.or", "mae.or")


save(list = c(beta_res, or_res), file = "sur-candes-p2.RData")
save.image(file = "full_sim-p2.RData")


