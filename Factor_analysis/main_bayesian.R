library(pheatmap)
source("Factor_analysis/Methods/bayesian_methods.R")

# -- load and process data --
load("Factor_analysis/Macro_data/stockwatson.Rda")
X <- stockwatson[-c(1:2, 197:200),] #n x p

# -- initialize --

nrun <- 50000
burnin <- 0
thin <- 100
k_tilde_0 <- floor(2*log(dim(X)[2]))
n <- dim(X)[1]
p <- dim(X)[2]
save_folder <- "Iterations"
# save_folder <- "Iterations_extended"
seeds <- c(1,2) 



# -- prob of adapting --
t <- seq(1:10000)
plot(t,1/exp(1+5e-4*t), type = 'l', ylab="Probability of Adapting")

# -- run gibbs sampler --
run_gibbs(X=X, nrun=nrun, thin=thin, k_tilde=k_tilde_0,
  save_folder=save_folder, seeds=seeds)
  
  
  # -- assess k_tilde and k_star --
k_func <- assess_ks(nrun=nrun, thin=thin, save_folder=save_folder)
k_star_post <- k_func$k_star_post
n_samples <- k_func$n_samples


# -- assess lambda multimodel post --
mcmc <- create_lambda_object(nrun=nrun, burnin=burnin, thin=thin,
                             k_star_post=k_star_post, n_samples=n_samples,
                             save_folder=save_folder, n=n, p=p, LQ=FALSE)
n_samples <- nrow(mcmc$L_mcmc_1)
assess_L_i(mcmc=mcmc, param_i=78, n_samples=n_samples, type="Lambda")


# -- assess diagnostics for L identifiable matrix --
n_samples <- k_func$n_samples
mcmc <- create_lambda_object(nrun=nrun, burnin=burnin, thin=thin,
                             k_star_post=k_star_post, n_samples=n_samples,
                             save_folder=save_folder, n=n, p=p, LQ=TRUE)

n_samples <- nrow(mcmc$L_mcmc_1)

#10 worst perfomring params
assess_L_i(mcmc=mcmc, param_i=657, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=662, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=698, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=616, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=534, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=575, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=493, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=739, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=497, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=411, n_samples=n_samples, type="L")
assess_L_i(mcmc=mcmc, param_i=11, n_samples=n_samples, type="L")


# -- assess sigma convergence --
sigma_mcmc <- assess_sigma(nrun=nrun, burnin=burnin, thin=thin, k_star_post=k_star_post,
                           n_samples=n_samples, save_folder=save_folder)
diagnostics(sigma_mcmc[,,c(17, 22, 21)])



# -- loadings interpretation --
macro_desc <- read_excel("Factor_analysis/Macro_data/macro_desc.xlsx")
L_post <- colMeans(rbind(mcmc$L_mcmc_1, mcmc$L_mcmc_2))
L_post <- matrix(L_post, nrow=p, ncol=k_star_post)
rownames(L_post) <- macro_desc$`Short name`
View(L_post)

pheatmap(L_post, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "L Matrix Heatmap",
         display_numbers = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)




