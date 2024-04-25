library(ggplot2)
set.seed(0)

# -- Generate data --
variable_1 <- rnorm(100)
variable_2 <- variable_1 * 0.5 + rnorm(100)
df <- data.frame(Variable1 = variable_1, Variable2 = variable_2)

# -- Perform PCA -- 
X <- scale(as.matrix(df))
SVD <- svd(X)
D <- SVD$d
theta <- SVD$v
pc1 <- theta[, 1] * 3
pc2 <- theta[, 2] * 2

lines_df <- data.frame(
  x = c(-pc1[1], -pc2[1]),
  y = c(-pc1[2], -pc2[2]),
  xend = c(pc1[1],pc2[1]),
  yend = c(pc1[2],pc2[2]),
  label = c("PC1", "PC2")
)

# -- Plot PCA and data points --
p <- ggplot(df, aes(x = Variable1, y = Variable2)) +
  geom_point(color = "blue") +
  theme_minimal()+
  coord_fixed()+
  xlab("Var 1") +  
  ylab("Var 2")

p <- p + geom_segment(data = lines_df, aes(x = x, y = y, xend = xend, yend = yend),
                      size = 1, color = "black")
          
p + geom_text(data = lines_df, aes(x = x, y = y, label = label), 
              vjust = "inward", hjust = "inward", nudge_x = 0.4, nudge_y = 0)
  



