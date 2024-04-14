library(readxl)
library(tensr)
library(coda)
source("Factor_analysis/Methods/graphical_diagnostics.R")
source("Factor_analysis/Methods/gibbs_sampler.R")

run_gibbs <- function(X, nrun, thin, k_tilde_0, save_folder, seeds){
  
  for(chain in 1:2){
    set.seed(seeds[chain])
    print(paste0("CHAIN ",chain))
    ks_t <- gibbs(data=X, k_tilde=k_tilde_0, nrun=nrun, chain=chain, thin=thin,
                  save_folder=save_folder)
    k_tilde_t <- ks_t$k_tilde_t
    k_star_t <- ks_t$k_star_t
    # -- Save k_tilde and k_star --
    file_name <- paste0("Factor_analysis/Results/",save_folder,"/chain_",chain,"/k/k_tilde.RData")
    save(k_tilde_t, file = file_name)
    file_name <- paste0("Factor_analysis/Results/",save_folder,"/chain_",chain,"/k/k_star.RData")
    save(k_star_t, file = file_name)
  }
}


assess_ks <- function(nrun, thin, save_folder){
  
  k_tilde_samples <- array(NA, c(nrun,2,1))
  dimnames(k_tilde_samples) <-  list(NULL, NULL, paste("k_tilde", sep=""))
  k_star_samples <- array(NA, c(nrun,2,1))
  dimnames(k_star_samples) <-  list(NULL, NULL, paste("k_star", sep=""))
  
  for(chain in 1:2){
      file_name <-paste0("Factor_analysis/Results/", save_folder,"/chain_",chain,"/k/k_tilde.RData")
      load(file = file_name)
      k_tilde_samples[,chain,] <- k_tilde_t
      
      file_name <-paste0("Factor_analysis/Results/", save_folder,"/chain_",chain,"/k/k_star.RData")
      load(file = file_name)
      k_star_samples[,chain,] <- k_star_t
    }
  #assess k_tilde
  plot(ts(k_tilde_samples[,1,]), ylab="k_tilde")
  lines(ts(k_tilde_samples[,2,]), col="red")
  
  #find optimal k_star
  tbl <- table(k_star_samples)
  k_star_post <- as.numeric(names(tbl)[which.max(tbl)])
  plot(ts(k_star_samples[,1,]), ylab="k_star")
  lines(ts(k_star_samples[,2,]), col="red")
  
  #assess convergence k_star
  print(sum(apply(k_star_samples, c(3), effectiveSize)))
  diagnostics(k_star_samples)
  
  #same samples as kept in Lambda after thinning
  k_tilde_samples <- k_tilde_samples[seq(from=thin+1, to=nrun, by=thin),,]
  n_samples <- min(colSums(k_star_post<=k_tilde_samples))
  return(list(k_star_post=k_star_post, n_samples=n_samples))
}

create_lambda_object <- function(nrun, burnin, thin, k_star_post,
                                 n_samples, save_folder, n, p){
  
  # -- select lambda with dim greater than k_star -- 
  L_mcmc_1 <- matrix(NA, nrow=n_samples, ncol=k_star_post*p)
  L_mcmc_2 <- matrix(NA, nrow=n_samples, ncol=k_star_post*p)
  i <- 1
  j <- 1
  for(t in (thin+1):(nrun/thin)){
    if((i>n_samples) & (j>n_samples)){
      break
    }
    iter <- format(t, scientific = FALSE)
    #chain 1
    if(i<=n_samples){
      file_name <-paste0("Factor_analysis/Results/",save_folder,"/chain_",1,"/Lambda/iter_",iter,".RData")
      load(file = file_name)
      if(ncol(lambda)>=k_star_post){
        # -- Apply LQ Decomp
        L <- lq(lambda)$L
        
        #truncate to effective number of factors
        L <- L[,1:k_star_post]
       
        L_mcmc_1[i,] <- as.array(L)
        i <- i+1
      }
    }
    
    #chain 2
    if(j<=n_samples){
      file_name <-paste0("Factor_analysis/Results/",save_folder,"/chain_",2,"/Lambda/iter_",iter,".RData")
      load(file = file_name)
      if(ncol(lambda)>=k_star_post){
        # -- Apply LQ Decomp
        L <- lq(lambda)$L
        
        #truncate to effective number of factors
        L <- L[,1:k_star_post]
        
        L_mcmc_2[j,] <- as.array(L)
        j <- j+1
      }
    }
    
  }
  #remove na rows
  L_mcmc_1 <- na.omit(L_mcmc_1)
  L_mcmc_2 <- na.omit(L_mcmc_2)
  n_samples <- min(nrow(L_mcmc_1), nrow(L_mcmc_2)) 
  L_mcmc_1 <- L_mcmc_1[1:n_samples,]  
  L_mcmc_2 <- L_mcmc_2[1:n_samples,]  
  
  
  
  #effective sample size
  ess <- effectiveSize(L_mcmc_1) + effectiveSize(L_mcmc_2)
  print(head(sort(ess[ess>0])))
  return(list(L_mcmc_1=L_mcmc_1, L_mcmc_2=L_mcmc_2))
}

assess_L_i <- function(mcmc, param_i, n_samples){

  #assess convergence of param_i
  L_i <- array(NA, c(n_samples,2,1))
  dimnames(L_i) <-  list(NULL, NULL, paste("L_", param_i, sep=""))
  L_i[,1,] <- mcmc$L_mcmc_1[1:n_samples,param_i]
  L_i[,2,] <- mcmc$L_mcmc_2[1:n_samples,param_i]
  print(sum(apply(L_i, c(3), effectiveSize)))
  diagnostics(L_i)
}

assess_sigma <- function(nrun, k_star_post, burnin, p, n_samples, save_folder){
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
      file_name <-paste0("Factor_analysis/Results/",save_folder,"/chain_",1,"/_iter_",n,".RData")
      load(file = file_name)
      if(k_star==k_star_post){
        sigma_mcmc[i,1,] <- as.array(sigma)
        i <- i+1
      }
    }
    
    #chain 2
    if(j<=n_samples){
      file_name <-paste0("Factor_analysis/Results/",save_folder,"/chain_",2,"/_iter_",n,".RData")
      load(file = file_name)
      if(k==k_star_post){
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