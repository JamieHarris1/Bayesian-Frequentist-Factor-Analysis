library(mvtnorm)
library(ggplot2)

set.seed(123)
# Number of observations
n <- 1000

# Number of observed variables
p <- 5

# Define factor loadings
loadings <- matrix(runif(2 * p), ncol = 2)

# Generate random latent factors
latent_factors <- rmvnorm(n, mean = c(0, 0), sigma = matrix(c(2, 0.8, 0.8, 1), nrow = 2, byrow = TRUE))

# F_latent (qxn)
F_latent <- t(latent_factors)

# Create observed variables
# X (pxn)
X <- loadings %*% F_latent + matrix(rnorm(p * n), nrow = p, ncol = n)

# Perform PCA
S <- cor(t(X))
eig <- eigen(S)
theta <- eig$vectors
alpha <- eig$values
variance <- alpha * 100 / sum(alpha)
cumvar <- cumsum(variance)
pca <- data.frame(
  eig = alpha,
  variance = variance,
  cumvariance = cumvar
)

# Create a scree plot
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


# Display the scree plot
print(scree_plot)
# ggsave("scree_plot.png", scree_plot, width = 6, height = 4)
