rm(list = ls())
set.seed(1991)
n <- 1000
p <- 200
b0 = c(3,1.5,0, -1.5,-3)
beta = rep(b0,each = p/5)
X <- matrix(rnorm(n * (p), 0, sqrt(1 / n)), n, p)
y <- rbinom(n, 1, plogis(X %*% beta))

m <- nrow(X)
y_dy <- p / (p + m) * 0.5 + m / (p + m) * y

library(microbenchmark)
library(RcppNumerical)
library(brglm2)
microbenchmark("DY" = fastLR(X, y_dy),
               "FIRTH" =  glm(y ~ -1 + X,
                              family = binomial("logit"),
                              method = "brglmFit", type = "AS_mean"),unit = "ms")
sessionInfo()




# M1 - openblas
# R version 4.1.1 (2021-08-10)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Big Sur 11.6.2
#
# Matrix products: default
# BLAS:   /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRblas.openblas.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/lib/libRlapack.dylib
# Unit: milliseconds
# expr        min         lq       mean     median         uq  max neval cld
# DY   1.132666   1.144433   1.191476   1.203145   1.221595  1.3735   100  a
# FIRTH 522.520892 531.041860 574.314200 536.368150 576.228985 877.7474   100   b


# M1 - standard BLAS library
# R version 4.1.1 (2021-08-10)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Big Sur 11.6.2
# Unit: milliseconds
# expr         min          lq        mean      median         uq         max neval cld
# DY    1.135413    1.148615    1.208988    1.209213    1.22346    2.757537    100  a
# FIRTH 1033.778879 1052.329616 1080.138944 1083.305444 1096.14328 1199.019047    100   b

## Consider larger values of N
rm(list = ls())
set.seed(1991)
n <- 10000
p <- 2000
b0 = c(3,1.5,0, -1.5,-3)
beta = rep(b0,each = p/5)
X <- matrix(rnorm(n * (p), 0, sqrt(1 / n)), n, p)
y <- rbinom(n, 1, plogis(X %*% beta))

m <- nrow(X)
y_dy <- p / (p + m) * 0.5 + m / (p + m) * y

library(RcppNumerical)
library(brglm2)
# Focus on a single replication
fastLR(X, y_dy)
# Takes a considerable amount of time
# brglm_fit(X,y, family = binomial("logit"), control = list(type = "AS_mean") )