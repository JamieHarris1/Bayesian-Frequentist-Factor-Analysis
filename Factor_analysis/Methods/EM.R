em <- function(X, k, iter){
  p <- dim(X)[1]
  n <- dim(X)[2]
  Sigma_tilde <- 1/n * X %*% t(X)
  
  # --- Initialise ---
  Lambda <- matrix(svd(t(X))$v[,1:k], nrow = p, ncol = k)
  Sigma <- diag(rep(1,p))
  I <- diag(rep(1,k))
  
  # --- EM Algorithm ---
  for(t in 1:iter){
    if(t%%1000==0){
      print(paste0("Iteration: ",t))
    }
    #E-STEP
    Sigma_inv <- solve(Sigma)
    M <- solve(I + t(Lambda) %*% Sigma_inv %*% Lambda)
    S <- M %*% t(Lambda) %*% Sigma_inv %*% X
    A <- n * M + S %*% t(S)
    
    #M-STEP
    Lambda_new <- X %*% t(S) %*% solve(A)
    Sigma_new <- diag(diag(Sigma_tilde - Sigma_tilde %*% Sigma_inv %*% Lambda %*% M %*% t(Lambda_new)))
    
    #UPDATE
    Sigma <- Sigma_new
    Lambda <- Lambda_new
  }
  
  return(list(Lambda=Lambda, Sigma=Sigma))
}