library(statmod)
library(matlab)
library(tensr)

gibbs <- function(data, k, nrun, chain, thin, save_folder)
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
  lambda <- matrix(0, nrow=p, ncol=k_initial)
  sigma <- rgamma(p, shape = asigma, scale = 1 / bsigma)
  Sigma <- diag(1 / sigma)
  eta <-  matrix(rnorm(n * k), n, k)
  mean_F <- zeros(n, k)
  var_F <- eye(k)
  phijh <- matrix(rgamma(p * k, shape = nu / 2, scale = 2 / nu), p, k)
  delta <- rgamma(k, shape=c(ad1, rep(ad2, k-1)), scale=c(bd1, rep(bd2, k-1)))
  tauh <- cumprod(delta)
  precision_lambda <- matvec(phijh, tauh)
  
  # -- Gibbs sampler -- 
  for(i in 1:nrun)
  {
    if(i%%1000==0){
      print(paste0("Iteration: ",i))
    }
    
    
    # -- update factors -- 
    sigma_lambda <- vecmat(sigma, lambda)
    var_F_inv <- diag(k) + t(sigma_lambda) %*% lambda
    T <- chol(var_F_inv)
    qrT <- qr(T)
    Q <- qr.Q(qrT)
    R <- qr.R(qrT)
    S <- solve(R)
    var_F <- tcrossprod(S)
    mean_F <- X %*% sigma_lambda %*% var_F
    x <- matrix(rnorm(n * k), nrow = n, ncol = k)
    F <- mean_F + x %*% t(S)
    
    # -- update lambda --
    F2 <- crossprod(F)
    for(j in 1:p)
    {
      Qlam <- diag(precision_lambda[j,]) + sigma[j] * F2
      blam <- sigma[j] * (t(F) %*% X[,j])
      Llam <- t(chol(Qlam))
      zlam <- rnorm(k)
      vlam <- forwardsolve(Llam, blam)
      mlam <- backsolve(t(Llam), vlam)
      ylam <- backsolve(t(Llam), zlam)
      lambda[j,] <- t(ylam + mlam)
    }
    
    # -- assess truncation --
    prob <- 1 / exp(b0 + b1 * i)
    lind <- colSums(as.matrix(abs(lambda) < epsilon)) / p
    vec <- lind >= prop
    num <- sum(vec)
    
    # -- update phi --
    for(h in 1:k)
      phijh[, h] <- rgamma(p, shape= (nu + 1) / 2, rate = (nu + tauh[h] * lambda[,h]^2) / 2)
    
    # --update tau --
    mat <- phijh * lambda^2
    ad <- ad1 + 0.5 * p * k
    bd <- bd1 + 0.5 * (1 / delta[1])  * sum(tauh * colSums(mat))
    delta[1] <- rgamma(1, shape = ad, scale = 1 / bd)
    tauh <- cumprod(delta)
    for(h in 2:k)
    {
      ad <- ad2 + 0.5 * p * (k - h + 1)
      bd <- bd2 + 0.5 * (1 / delta[h]) * sum(tauh[h:k] * colSums(as.matrix(mat[, h:k])))
      delta[h] <- rgamma(1, shape = ad, scale = 1 / bd)
      tauh <- cumprod(delta)
    }
    
    # -- update sigma --
    Xtil <- X - F %*% t(lambda)
    sigma <- rgamma(p, shape = asigma + 0.5 * n, rate = bsigma + 0.5 * colSums(Xtil^2))
    Sigma <- diag(1 / sigma)
    
    # -- update lambda precision
    precision_lambda <- matvec(phijh, tauh)
    
    #Apply thinning and save
    if(i %% thin == 0) {
      file_name <- paste0("Factor_analysis/Results/",save_folder,"/chain_",chain,"_iter_",i,".RData")
      save(k, sigma, lambda , file = file_name)
    }
    
    
    # --- Make adaptations ---
    uu <- runif(1, min=0, max=1)
    if (uu < prob) {
      if (i > 0 && num == 0 && all(lind < 0.995)) {
        k <-  k + 1
        lambda <- cbind(lambda, rep(0, p))
        F <- cbind(F, rnorm(n, 0, 1))
        phijh <- cbind(phijh, rgamma(p, shape = nu/2, scale = 2/nu))
        delta[k] <- rgamma(1, shape = ad2, scale = 1/bd2)
        tauh <- cumprod(delta)
        precision_lambda <- matvec(phijh, tauh)
      }
      else if (num > 0) {
        # print("REMOVING FACTOR")
        nonred <- setdiff(1:k, which(vec))
        k <- max(k - num, 1)
        lambda <- lambda[, nonred]
        phijh <- phijh[, nonred]
        F <- F[, nonred]
        delta <- delta[nonred]
        tauh <- cumprod(delta)
        precision_lambda <- matvec(phijh, tauh)
      }
    }
    
    
    
    
  }
}