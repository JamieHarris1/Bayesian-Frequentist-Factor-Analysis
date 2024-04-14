source("Factor_analysis/Methods/bayesian_methods.R")

# -- load and process data --
load("Factor_analysis/Macro_data/stockwatson.Rda")
X <- stockwatson[-c(1:2, 197:200),] #n x p

# -- initialize --
nrun <- 3000000
burnin <- 0
thin <- 100
k_tilde_0 <- floor(2*log(dim(X)[2]))
n <- dim(X)[1]
p <- dim(X)[2]
save_folder <- "Iterations"
seeds <- c(1,2) 


# -- prob of adapting --
t <- seq(1:10000)
plot(t,1/exp(1+5e-4*t), type = 'l', ylab="Probability of Adapting")

# -- run gibbs sampler --
# run_gibbs(X=X, nrun=nrun, thin=thin, k_tilde=k_tilde_0,
#           save_folder=save_folder, seeds=seeds)


# -- assess k_tilde and k_star --
k_func <- assess_ks(nrun=nrun, thin=thin, save_folder=save_folder)
k_star_post <- k_func$k_star_post
n_samples <- k_func$n_samples



# -- assess lambda convergence --
mcmc <- create_lambda_object(nrun=nrun, burnin=burnin, thin=thin, k_star_post=k_star_post,
                             n_samples=n_samples, save_folder=save_folder, n=n, p=p)

n_samples <- nrow(mcmc$L_mcmc_1)
assess_L_i(mcmc=mcmc, param_i=552, n_samples=n_samples)
assess_L_i(mcmc=mcmc, param_i=547, n_samples=n_samples)
assess_L_i(mcmc=mcmc, param_i=657, n_samples=n_samples)
assess_L_i(mcmc=mcmc, param_i=662, n_samples=n_samples)
assess_L_i(mcmc=mcmc, param_i=616, n_samples=n_samples)

# -- loadings interpretation --
macro_desc <- read_excel("Factor_analysis/Macro_data/macro_desc.xlsx")
lambda_post <- colMeans(rbind(mcmc$lambda_mcmc_1, mcmc$lambda_mcmc_2))
lambda_post <- matrix(lambda_post, nrow=p, ncol=k_star)
rownames(lambda_post) <- macro_desc$`Short name`
View(lambda_post)

# -- assess sigma convergence --
assess_sigma(nrun=nrun, k_star_post=k_star_post, burnin=burnin, p=p, n_samples=n_samples, save_folder=save_folder)


