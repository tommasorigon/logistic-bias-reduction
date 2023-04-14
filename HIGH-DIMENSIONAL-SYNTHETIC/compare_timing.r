rm(list = ls())
set.seed(1991)
n = 1000
p = 200
b0 = c(3,1.5,0, -1.5,-3)
beta = rep(b0,each = p/5)
X = matrix(rnorm(n * (p), 0, sqrt(1 / n)), n, p)
y = rbinom(n, 1, plogis(X %*% beta))

m = nrow(X)
y_dy = p / (p + m) * 0.5 + m / (p + m) * y

library(microbenchmark)
library(RcppNumerical)
library(brglm2)
times1 <- microbenchmark("DY" = fastLR(X, y_dy),
               "ML" = fastLR(X, y),
               "DY_glm" = glm(y_dy ~ -1 + X, family = binomial("logit")),
               "ML_glm" = glm(y ~ -1 + X, family = binomial("logit")),
               "FIRTH" =  glm(y ~ -1 + X, family = binomial("logit"),  method = "brglmFit", type = "AS_mean"),
               unit = "ms")

times1

#Unit: milliseconds
#expr         min          lq        mean      median          uq         max neval
#DY    1.417825    1.536040    1.657501    1.618864    1.734205    2.675760   100
#ML    1.515803    1.637326    1.746349    1.702804    1.828086    2.472775   100
#DY_glm  111.607941  119.617553  126.308349  123.216071  127.418655  259.807011   100
#ML_glm  138.646314  145.807803  150.875042  149.376876  153.386644  233.923078   100
#FIRTH 1045.248981 1085.919382 1168.210768 1183.048705 1204.458069 1479.942433   100

sessionInfo()


## Consider a settings with correlated covariates. Firth takes more than 2 hours
n = 10000
p = 2000
b0 = c(3,1.5,0, -1.5,-3)
beta = rep(b0,each = p/5)
S = tcrossprod(matrix(rnorm(p*10),p)) + diag(1/n,p,p)
X = mvtnorm::rmvnorm(n, sigma = S)
y = rbinom(n, 1, plogis(X %*% beta))
m = nrow(X)
y_dy = p / (p + m) * 0.5 + m / (p + m) * y


# Focus on a single replication (takes roughly 1 second)
fastLR(X, y_dy)
# Takes a considerable amount of time (> 2 hours)
t0 = Sys.time()
brglm_fit(X,y, family = binomial("logit"), control = list(type = "AS_mean",maxit = 5e3))
t1 = t0 - Sys.time()
