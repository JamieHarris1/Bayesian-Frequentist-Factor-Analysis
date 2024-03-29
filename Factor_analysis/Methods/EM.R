em <- function(data, k, iter){
  x <- as.matrix(data)
  x <- matrix(as.numeric(x), nrow = nrow(x)) #nxp
  x <- t(x) #pxn
  p <- dim(x)[1]
  n <- dim(x)[2]
  cov_hat <- 1/n * x %*% t(x)
  epsilon <- 1e-02
  
  #PCA lam and random psi
  lambda <- matrix(eigen(cov_hat)$vectors[,1:k], nrow = p, ncol = k)
  sigma <- diag(runif(p))
  I <- diag(rep(1,k))
  
  #loglikelihood
  c <- (n*p/2) * log(2*pi)
  LL <- rep(0, 1)
  
  #####
  
  for(j in 1:iter){
    #E-STEP
    sigma_inv <- solve(sigma)
    M <- solve(I + t(lambda) %*% sigma_inv %*% lambda)
    S <- M %*% t(lambda) %*% sigma_inv %*% x
    A <- n * M + S %*% t(S)
    
    #M-STEP
    lambda_new <- x %*% t(S) %*% solve(A)
    sigma_new <- diag(diag(cov_hat - cov_hat %*% sigma_inv %*% lambda %*% M %*% t(lambda_new)))
    
    #UPDATE
    sigma <- sigma_new
    lambda <- lambda_new
  }
  
  return(list(lambda=lambda, sigma=sigma))
}
