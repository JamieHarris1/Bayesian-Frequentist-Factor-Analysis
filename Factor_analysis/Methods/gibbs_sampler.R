library(statmod)
library(matlab)
library(tensr)

gibbs <- function(data, k_tilde, nrun, chain, thin, save_folder)
{
  # -- read data --
  p <- ncol(data)
  n <- nrow(data)
  X <- scale(data)
  
  # -- hyper-parameters --
  asigma <- 1
  bsigma <- 0.3
  nu <- 3
  epsilon <- 1e-3
  prop <- 1.00
  b0 <- 1
  b1 <- 0.0001
  ad1 <- 2.1
  ad2 <- 3.2
  bd1 <- 1
  bd2 <- 1
  
  # -- initial values --
  lambda <- matrix(0, nrow=p, ncol=k_tilde)
  sigma <- rgamma(p, shape = asigma, scale = 1 / bsigma)
  Sigma <- diag(1 / sigma)
  eta <-  matrix(rnorm(n * k_tilde), n, k_tilde)
  mean_F <- zeros(n, k_tilde)
  var_F <- eye(k_tilde)
  phijh <- matrix(rgamma(p * k_tilde, shape = nu / 2, scale = 2 / nu), p, k_tilde)
  delta <- rgamma(k_tilde, shape=c(ad1, rep(ad2, k_tilde-1)), scale=c(bd1, rep(bd2, k_tilde-1)))
  tauh <- cumprod(delta)
  precision_lambda <- matvec(phijh, tauh)
  k_star <- k_tilde_0
  
  # -- Set up output --
  k_tilde_t <- c()
  k_star_t <- c()
  
  
  # -- Gibbs sampler -- 
  for(t in 1:nrun)
  {
    if(t%%1000==0){
      print(paste0("Iteration: ",t))
    }
    
    
    # -- update factors -- 
    sigma_lambda <- vecmat(sigma, lambda)
    var_F_inv <- diag(k_tilde) + t(sigma_lambda) %*% lambda
    T <- chol(var_F_inv)
    qrT <- qr(T)
    Q <- qr.Q(qrT)
    R <- qr.R(qrT)
    S <- solve(R)
    var_F <- tcrossprod(S)
    mean_F <- X %*% sigma_lambda %*% var_F
    x <- matrix(rnorm(n * k_tilde), nrow = n, ncol = k_tilde)
    F <- mean_F + x %*% t(S)
    
    # -- update lambda --
    F2 <- crossprod(F)
    for(j in 1:p)
    {
      Qlam <- diag(precision_lambda[j,]) + sigma[j] * F2
      blam <- sigma[j] * (t(F) %*% X[,j])
      Llam <- t(chol(Qlam))
      zlam <- rnorm(k_tilde)
      vlam <- forwardsolve(Llam, blam)
      mlam <- backsolve(t(Llam), vlam)
      ylam <- backsolve(t(Llam), zlam)
      lambda[j,] <- t(ylam + mlam)
    }
    
    # -- update phi --
    for(h in 1:k_tilde)
      phijh[, h] <- rgamma(p, shape= (nu + 1) / 2, rate = (nu + tauh[h] * lambda[,h]^2) / 2)
    
    # --update tau --
    mat <- phijh * lambda^2
    ad <- ad1 + 0.5 * p * k_tilde
    bd <- bd1 + 0.5 * (1 / delta[1])  * sum(tauh * colSums(mat))
    delta[1] <- rgamma(1, shape = ad, scale = 1 / bd)
    tauh <- cumprod(delta)
    for(h in 2:k_tilde)
    {
      ad <- ad2 + 0.5 * p * (k_tilde - h + 1)
      bd <- bd2 + 0.5 * (1 / delta[h]) * sum(tauh[h:k_tilde] * colSums(as.matrix(mat[, h:k_tilde])))
      delta[h] <- rgamma(1, shape = ad, scale = 1 / bd)
      tauh <- cumprod(delta)
    }
    
    # -- update sigma --
    Xtil <- X - F %*% t(lambda)
    sigma <- rgamma(p, shape = asigma + 0.5 * n, rate = bsigma + 0.5 * colSums(Xtil^2))
    Sigma <- diag(1 / sigma)
    
    # -- update lambda precision
    precision_lambda <- matvec(phijh, tauh)
    
    
    #Save every k_tilde and k_star
    k_tilde_t[t] <- k_tilde
    k_star_t[t] <- k_star
    
    #Save thinned Lambda and Sigma
    if(t %% thin == 0) {
      lambda_file <- paste0("Factor_analysis/Results/",save_folder,"/chain_",chain,"/Lambda/iter_",(t/thin),".RData")
      save(lambda , file = lambda_file)
      sigma_file <- paste0("Factor_analysis/Results/",save_folder,"/chain_",chain,"/sigma/iter_",(t/thin),".RData")
      save(sigma , file = sigma_file)
    }
    
    # -- number of effective factors --
    lind <- colSums(as.matrix(abs(lambda) < epsilon)) / p
    vec <- lind >= prop
    m_t <- sum(vec)
    k_star <- k_tilde - m_t
    
    # --- Make adaptations ---
    prob <- 1 / exp(b0 + b1 * t)
    uu <- runif(1, min=0, max=1)
    if (uu < prob) {
      if (t > 0 && m_t == 0 && all(lind < 0.995)) {
        k_tilde <-  k_tilde + 1
        lambda <- cbind(lambda, rep(0, p))
        F <- cbind(F, rnorm(n, 0, 1))
        phijh <- cbind(phijh, rgamma(p, shape = nu/2, scale = 2/nu))
        delta[k_tilde] <- rgamma(1, shape = ad2, scale = 1/bd2)
        tauh <- cumprod(delta)
        precision_lambda <- matvec(phijh, tauh)
      }
      else if (m_t > 0) {
        # print("REMOVING FACTOR")
        nonred <- setdiff(1:k_tilde, which(vec))
        k_tilde <- max(k_star, 1)
        lambda <- lambda[, nonred]
        phijh <- phijh[, nonred]
        F <- F[, nonred]
        delta <- delta[nonred]
        tauh <- cumprod(delta)
        precision_lambda <- matvec(phijh, tauh)
      }
    }
    
    
    
    
  }
  return(list(k_tilde_t=k_tilde_t, k_star_t=k_star_t))
}