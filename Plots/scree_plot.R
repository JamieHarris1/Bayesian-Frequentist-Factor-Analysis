library(mvtnorm)
library(ggplot2)

set.seed(123)

# --- Generate Data ---
p <- 5
n <- 1000
loadings <- matrix(runif(2 * p), ncol = 2)
latent_factors <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(2, 0.8, 0.8, 1), nrow = 2, byrow = TRUE))
X <- loadings %*% t(latent_factors) + matrix(rnorm(p * n), nrow = p, ncol = n)

# --- Perform PCA ---
Sx <- cov(t(X))
eig <- eigen(Sx)
theta <- eig$vectors
alpha <- eig$values

# SVD <- svd(t(X))
# D <- SVD$d
# alpha <- D^2 / (n-1)

variance <- alpha * 100 / sum(alpha)
cumvar <- cumsum(variance)
pca <- data.frame(
  eig = alpha,
  variance = variance,
  cumvariance = cumvar
)




# --- Create Scree Plot ---
scree_plot <- ggplot(pca, aes(x = 1:nrow(pca), y = variance)) +
  geom_line(stat = "identity", color = "skyblue", size = 1.5) +
  geom_point(aes(x = 1:nrow(pca), y = variance), color = "red", size = 1) +
  geom_text(aes(x = 1:nrow(pca), y = variance, label = round(variance, 2)),
            vjust = -0.5, hjust = 0, color = "black", size = 3) +
  labs(title = "Variance explained by each PC",
       x = "Principal Component - j",
       y = expression(V[j]),
       )+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA),  
        panel.background = element_rect(fill = "white"),
        panel.grid = element_blank(),
        panel.border = element_blank() )

print(scree_plot)

