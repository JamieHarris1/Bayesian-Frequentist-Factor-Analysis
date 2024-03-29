library(readxl)
library(tensr)
library(coda)
source("Factor_analysis/Methods/graphical_diagnostics.R")
source("Factor_analysis/Methods/gibbs_sampler.R")

run_gibbs <- function(X, nrun, thin, k_initial, save_folder, seed1, seed2){
  set.seed(seed1)
  gibbs(data=X, k=k_initial, nrun=nrun, chain=1, thin=thin,
        save_folder=save_folder)
  
  print("CHAIN 2")
  
  set.seed(seed2)
  gibbs(data=X, k=k_initial, nrun=nrun, chain=2, thin=thin,
        save_folder=save_folder)
  
  
}

find_kstar <- function(nrun, thin, save_folder){
  # --- Find optimal num of factors: k_star ---\
  k_t <- matrix(NA, nrow=2, ncol=(nrun/thin))
  for(chain in 1:2){
    for(n in seq(from=thin, to=nrun, by=thin)){
      iter <- format(n, scientific = FALSE)
      file_name <-paste0("Factor_analysis/Results/", save_folder,"/chain_",chain,"_iter_",
                         iter,".RData")
      load(file = file_name)
      k_t[chain,n/thin] <- k
      
    }
  }
  tbl <- table(k_t)
  k_star <- as.numeric(names(tbl)[which.max(tbl)])
  plot(ts(k_t[1,]), ylab="k")
  lines(ts(k_t[2,]), col="red")
  return(list(k_star=k_star, k_t=k_t))
}

create_lambda_object <- function(nrun, thin, k_star, k_t, save_folder, n, p){
  # -- select lambda with dim k_star -- 
  n_samples <- min(rowSums(k_t==k_star))
  lambda_mcmc_1 <- matrix(NA, nrow=n_samples, ncol=k_star*p)
  lambda_mcmc_2 <- matrix(NA, nrow=n_samples, ncol=k_star*p)
  i <- 1
  j <- 1
  for(n in seq(from=thin, to=nrun, by=thin)){
    if((i>n_samples) & (j>n_samples)){
      break
    }
    iter <- format(n, scientific = FALSE)
    #chain 1
    if(i<=n_samples){
      file_name <-paste0("Factor_analysis/Results/",save_folder,"/chain_",1,"_iter_",iter,".RData")
      load(file = file_name)
      if(k==k_star){
        # -- Apply LQ Decomp
        lambda <- lq(lambda)$L
        lambda_mcmc_1[i,] <- as.array(lambda)
        i <- i+1
      }
    }
    
    #chain 2
    if(j<=n_samples){
      file_name <-paste0("Factor_analysis/Results/",save_folder,"/chain_",2,"_iter_",iter,".RData")
      load(file = file_name)
      if(k==k_star){
        # -- Apply LQ Decomp
        lambda <- lq(lambda)$L
        lambda_mcmc_2[j,] <- as.array(lambda)
        j <- j+1
      }
    }
    
  }
  
  #post thinning
  lambda_mcmc_1 <- lambda_mcmc_1[seq(from=1,to=n_samples,by=5),]
  lambda_mcmc_2 <- lambda_mcmc_2[seq(from=1,to=n_samples,by=5),]
  
  #effective sample size
  ess <- effectiveSize(lambda_mcmc_1) + effectiveSize(lambda_mcmc_2)
  print(head(sort(ess[ess>0])))
  return(list(lambda_mcmc_1=lambda_mcmc_1, lambda_mcmc_2=lambda_mcmc_2))
}

assess_lambda_i <- function(mcmc, param_i, n_samples){

  #assess convergence of param_i
  lambda_i <- array(NA, c(n_samples,2,1))
  dimnames(lambda_i) <-  list(NULL, NULL, paste("lambda_", param_i, sep=""))
  lambda_i[,1,] <- mcmc$lambda_mcmc_1[1:n_samples,param_i]
  lambda_i[,2,] <- mcmc$lambda_mcmc_2[1:n_samples,param_i]
  print(sum(apply(lambda_i, c(3), effectiveSize)))
  diagnostics(lambda_i)
}

assess_sigma <- function(nrun, burnin, p, n_samples, save_folder){
  #Sigma convergence
  i <- 1
  j <- 1
  sigma_mcmc <- array(NA, c(n_samples,2,p))
  for(n in seq(from=thin, to=nrun, by=thin)){
    if((i>n_samples) & (j>n_samples)){
      break
    }
    #chain 1
    if(i<=n_samples){
      file_name <-paste0("Factor_analysis/Results/",save_folder,"/chain_",1,"_iter_",n,".RData")
      load(file = file_name)
      if(k==k_star){
        sigma_mcmc[i,1,] <- as.array(sigma)
        i <- i+1
      }
    }
    
    #chain 2
    if(j<=n_samples){
      file_name <-paste0("Factor_analysis/Results/",save_folder,"/chain_",2,"_iter_",n,".RData")
      load(file = file_name)
      if(k==k_star){
        sigma_mcmc[j,2,] <- as.array(sigma)
        j <- j+1
      }
    }
    
  }
  dimnames(sigma_mcmc) <-  list(NULL, NULL, paste("sigma_", 1:p, sep=""))
  diagnostics(sigma_mcmc)
  ess <- sum(apply(sigma_mcmc, c(3), effectiveSize))
  head(sort(ess[ess>0]))
}