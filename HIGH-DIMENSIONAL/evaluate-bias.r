set.seed(1991)
n <- 500
p <- 100

SS <- matrix(0.2, p - 1, p - 1)
diag(SS) <- 1

X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1) %*% chol(SS))
betas <- runif(p)/10
y <- rbinom(n, 1, prob = plogis(X %*% betas))

m0 <- glm(y ~ X - 1, family = binomial("logit"))
Nsim = 1e4
y_sim = simulate(m0,nsim = Nsim)

ml <- br <- mbr <- clogg <- dy <- matrix(NA, nrow = Nsim, ncol = length(betas))

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

  # KENNE PAGUI ET AL. (2017)
  fit_kp <- glm(y ~ -1 + X,
                family = binomial("logit"),
                method = "brglmFit", type = "AS_median"
  )

  ml[i,] = coef(mle)
  dy[i,] = fit_fast_dy$coefficients
  clogg[i,] = fit_fast_clogg$coefficients
  br[i,] = coef(fit_firth)
  mbr[i,] = coef(fit_kp)
}
