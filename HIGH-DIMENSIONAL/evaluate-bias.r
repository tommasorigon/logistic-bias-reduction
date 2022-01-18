set.seed(1991)
n <- 1000
p <- 200
beta <- c(rep(10, p/8), rep(-10, p/8), rep(0, 3*p/4))
X <- matrix(rnorm(n * p, 0, sqrt(1/n)), n, p)
y <- rbinom(n, 1, plogis(X %*% beta))




m0 <- glm(y ~ X - 1, family = binomial("logit"))
Nsim = 1e3
y_sim = simulate(m0,nsim = Nsim)

ml <- br <- clogg <- dy <- matrix(NA, nrow = Nsim, ncol = length(beta))

library(brglm2)
library(RcppNumerical)
m <- nrow(X)
for(i in 1:Nsim) {

  if (i %% 10 == 0) print(i)
  y = y_sim[,i]

  mle <- glm(y ~ X - 1, family = binomial("logit"))

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



  ml[i,] = coef(mle)
  dy[i,] = fit_fast_dy$coefficients
  clogg[i,] = fit_fast_clogg$coefficients
  br[i,] = coef(fit_firth)
}


## BIAS in beta parameterization
bias.beta <- data.frame(
  ml = colMeans(ml[ ]),
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
bias.beta
df = reshape2::melt(bias.beta)
head(df)
library(ggplot2)
table(df$variable)
ggplot(df) +
  geom_boxplot(aes(y = value, color = variable, group = variable)) +
  theme_bw()

matplot(bias.beta)
matplot(mse.beta)
