## Contour plots

``` r
rm(list = ls())
set.seed(1991)
x1 <- rnorm(10)
x2 <- rnorm(10)
X <- cbind(x1, x2)
X <- scale(X)

dy_prior <- function(beta1, beta2, X, kappai = rep(0.5, NROW(X)), tau = NCOL(X) / NROW(X)) {
  nu <- X %*% c(beta1, beta2)
  out <- tau * (t(kappai) %*% nu - sum(log(1 + exp(nu))))
  return(exp(out))
}

jeff_prior <- function(beta1, beta2, X) {
  pis <- plogis(X %*% c(beta1, beta2))
  w <- pis * (1 - pis)
  XtX<- t(X) %*% diag(c(w), NROW(X), NROW(X)) %*% X
  return(exp(0.5* determinant(XtX, logarithm = TRUE)$modulus))
}

dgelman <- function(beta1, beta2) {
  dcauchy(beta1, 0, scale = 2.5) * dcauchy(beta2, 0, 2.5)
}

rr <- seq(-6, 6, l = 200)
gr <- gr2 <- gr3 <- expand.grid(rr, rr)

gr$values <- apply(gr, 1, function(x) dy_prior(x[1], x[2], X = X))
gr2$values <- apply(gr2, 1, function(x) jeff_prior(x[1], x[2], X = X))
gr3$values <- apply(gr3, 1, function(x) dgelman(x[1], x[2]))

# Normalizzazione stupida
library(cubature)

C_dy <- hcubature(function(x) dy_prior(x[1], x[2], X = X),c(-Inf,-Inf), c(Inf,Inf))$integral
C_jeff <- hcubature(function(x) jeff_prior(x[1], x[2], X = X),rep(-Inf,2), rep(Inf, 2))$integral

# Sanity check
# hcubature(function(x) dgelman(x[1], x[2]),rep(-Inf,2), rep(Inf, 2))$integral

gr$values <- gr$values / C_dy
gr2$values <- gr2$values / C_jeff

gr$w <- "Diaconis & Ylvisaker"
gr2$w <- "Jeffrey (Firth, 1993)"
gr3$w <- "Cauchy (Gelman et al., 2008)"
gr_c <- rbind(gr, gr2, gr3)

library(ggplot2)
library(metR)
p1 <- ggplot(gr_c) +
  geom_contour(aes(x = Var1, y = Var2, z = values), bins = 15, col = "gray") +
  geom_text_contour(aes(x = Var1, y = Var2, z = values), size = 2, rotate = TRUE) +
  facet_wrap(~w) +
  theme_bw() + xlab(expression(beta[1])) + ylab(expression(beta[2]))
p1
```

![](/Users/tommaso/Google%20Drive/University/Paper/R%20packages/logistic-bias-reduction/GRAPHS/README_files/figure-gfm/contour-1.png)<!-- -->

``` r
p2 <- ggplot(gr_c) +
  geom_contour_filled(aes(x = Var1, y = Var2, z = values), bins = 25) +
  facet_wrap(~w) +
  theme_bw() + xlab(expression(beta[1])) + ylab(expression(beta[2])) + theme(legend.position = "none") 
p2
```

![](/Users/tommaso/Google%20Drive/University/Paper/R%20packages/logistic-bias-reduction/GRAPHS/README_files/figure-gfm/contour-2.png)<!-- -->

``` r
plot(cars)
```

![](/Users/tommaso/Google%20Drive/University/Paper/R%20packages/logistic-bias-reduction/GRAPHS/README_files/figure-gfm/contour-3.png)<!-- -->
