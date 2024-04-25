library("tidyverse")

PCA <- function(X){
  # --- Perform PCA ---
  n <- dim(X)[2]
  print(n)
  SVD <- svd(t(X))
  Theta <- SVD$v
  print(dim(Theta))
  D <- SVD$d
  alpha <- D^2 / (n-1)
  
  variance <- alpha * 100 / sum(alpha)
  cumvar <- cumsum(variance)
  pca <- data.frame(
    eig = alpha,
    variance = variance,
    cumvariance = cumvar
  )
  
  # --- Create a Scree Plot ---
  scree_plot <- ggplot(pca, aes(x = 1:nrow(pca), y = variance)) +
    geom_line(stat = "identity", color = "skyblue", size = 1.5) +
    geom_point(aes(x = 1:nrow(pca), y = variance), color = "red", size = 1) +
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
  return(list(Theta=Theta, alpha=alpha))
}