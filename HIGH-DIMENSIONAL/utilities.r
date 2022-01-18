logit_NR <- function(X, y,  tol = 1e-16, beta_start = NULL, maxiter=10000){

  if (is.null(n <- nrow(X)))
    stop("'X' must be a matrix")

  loglik <- numeric(maxiter)

  # Initialization (If null, implicitely initialized at beta=0)
  if(is.null(beta_start)) {beta <- solve(crossprod(X/4,X),crossprod(X,y-0.5))} else {beta <- beta_start}
  # Initialization
  eta       <- c(X%*%beta)
  prob      <- 1/(1+exp(-eta))
  w         <- prob*(1-prob)

  # First value of the likelihood
  loglik[1] <- sum(y*eta - log(1+exp(eta)))

  # Iterative procedure
  for(t in 2:maxiter){
    beta       <- beta + solve(qr(crossprod(X*w,X)),crossprod(X,y-prob))
    eta        <- c(X%*%beta)
    prob       <- 1/(1+exp(-eta))
    w          <- prob*(1-prob)
    loglik[t]  <- sum(y*eta - log(1+exp(eta)))
    if(loglik[t] - loglik[t-1] < tol) return(list(beta=beta, Convergence=cbind(Iteration=(1:t)-1,
                                                                               Loglikelihood=loglik[1:t])))
  }
  stop("The algorithm has not reached convergence")
}