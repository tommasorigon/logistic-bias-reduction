set.seed(1)
n <- 1000
p <- 100

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

# H20 estimate
library(h2o)
h2o.init()

dataset <- as.h2o(data.frame(y = y, X = X[,-1]))
fit_h2o <- h2o.glm(family = "binomial",
                        x = colnames(dataset[,-1]),
                        y = "y",
                        training_frame = dataset,
                        lambda = 0,
                        compute_p_values = TRUE)


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

# # KENNE PAGUI ET AL. (2017)
# fit_kp <- glm(y ~ -1 + X,
#   family = binomial("logit"),
#   method = "brglmFit", type = "AS_median"
# )
# beta_kp <- coef(fit_kp)
# sd_kp <- summary(fit_kp)$coefficients[, 2]

c(elapsed_dy, elapsed_firth)
