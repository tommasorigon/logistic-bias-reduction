n_p <- lapply(as.list(c(1, 5, 10)), function(x) list(n = 100 * (10) * x, p = 10 * 10 * x))
disp_string <- "Starting with n=%s and p=%s\n"
library(brglm2)
library(fastglm)

out <- n_p
for (i in 1:length(n_p)) {
  set.seed(1)
  n <- n_p[[i]]$n
  p <- n_p[[i]]$p
  cat(sprintf(disp_string, n, p))

  SS <- matrix(0.2, p, p)
  diag(SS) <- 1

  X <- matrix(rnorm(n * p), n, p) %*% chol(SS)
  betas <- runif(p)
  y <- rbinom(n, 1, prob = plogis(X %*% betas))

  m <- nrow(X)
  y_clogg <- p / (p + m) * mean(y) + m / (p + m) * y
  y_dy <- p / (p + m) * 0.5 + m / (p + m) * y

  # FIRTH (1993)
  t0 <- Sys.time()
  fit_firth <- glm(y ~ X,
    family = binomial("logit"),
    method = "brglmFit", type = "AS_mean", control = brglmControl(maxit = 500, trace = T)
  )
  t1 <- Sys.time()
  elapsed_firth <- t1 - t0

  # fit_dy <- fastglm(x = X,y=y_dy,
  #  family = binomial("logit")
  # )
  XI <- model.matrix(fit_firth)
  t0 <- Sys.time()
  fit_dy2 <- fastglm(X, y, family = binomial())
  elapsed_dy2 <- Sys.time() - t0

  out[[i]]$time <- (list("dy" = elapsed_dy2, "firth" = elapsed_firth))
  out[[i]]$coef <- list(coef(fit_firth), coef(fit_dy2))

  cat(i, "\n")
  rm(list = c("X", "y", "fit_firth", "fit_dy2"))
  gc()
  save(out, file = "out.RData")
}
