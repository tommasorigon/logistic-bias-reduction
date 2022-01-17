set.seed(1)
n <- 2000
p <- 400

SS <- matrix(0.2, p - 1, p - 1)
diag(SS) <- 1

X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1) %*% chol(SS))
betas <- runif(p)
y <- rbinom(n, 1, prob = plogis(X %*% betas))

library(brglm2)
m <- nrow(X)
y_clogg <- p / (p + m) * mean(y) + m / (p + m) * y
y_dy <- p / (p + m) * 0.5 + m / (p + m) * y

# DY ESTIMATE
t0 <- Sys.time()
fit_dy <- glm(y_dy ~ X - 1, family = binomial("logit"))
t1 <- Sys.time()
elapsed_dy <- t1 - t0
beta_dy <- coef(fit_dy)
sd_dy <- summary(fit_dy)$coefficients[, 2]

library(RcppNumerical)
t0 <- Sys.time()
fit_fast_dy <- fastLR(X, y_dy)
t1 <- Sys.time()
elapsed_fast_dy <- t1 - t0
beta_fast_dy <- fit_fast_dy$coefficients


# FIRTH (1993)
t0 <- Sys.time()
fit_firth <- glm(y ~ -1 + X,
                 family = binomial("logit"),
                 method = "brglmFit", type = "AS_mean", control = brglmControl(maxit = 500, trace = T)
)
t1 <- Sys.time()
beta_firth <- coef(fit_firth)
sd_firth <- summary(fit_firth)$coefficients[, 2]
elapsed_firth <- t1 - t0

# KENNE PAGUI ET AL. (2017)
t0 <- Sys.time()
fit_kp <- glm(y ~ -1 + X,
  family = binomial("logit"),
  method = "brglmFit", type = "AS_median"
)
t1 <- Sys.time()
beta_kp <- coef(fit_kp)
elapsed_kp <- t1 - t0


c(elapsed_dy, elapsed_fast_dy, elapsed_firth, elapsed_kp)
