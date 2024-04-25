library(tidyverse)
library(pixmap)

plot_face <- function(mat){
  plot(pixmapGrey(255* matrix(mat, nrow=112, ncol=92)))
}

d_ij <- function(j,X,X_test, k){
  X_tilde_test <- t(Theta[,1:k]) %*% X_test[,j]
  d_i <- c()
  for(i in 1:dim(X)[2]){
    X_tilde <- t(Theta[,1:k]) %*% X[,i]
    Sx <- cov(X_tilde,X_tilde_test)
    d_i[i] <- sqrt(1/Sx * t(X_tilde-X_tilde_test) %*% (X_tilde-X_tilde_test))
  }
  id <- which.min(d_i)
  MH <- d_i[id]
  plot_face(X[,id])
  plot_face(X_test[,j])
  return(list(id=id, MH=MH))
}


options(warn = -1)
k <- 200
r <- matrix(NA, nrow=10304)

# --- Load Training Images ---
for(n in 11:409){
  if(!n%%10==0){
    file <- paste0("Facial_recognition/faces/s",n,".pgm")
    pgm_image <- read.pnm(file)
    I_i <- pgm_image@grey
    r <- cbind(r,c(I_i))
  }
}
r <- r[,-1]


# --- Compute Mean Face ---
Omega <- rowMeans(r)
X <- r - Omega

# --- Load Test Images ---
test <- matrix(NA, nrow=10304)
for(n in seq(110,4010,100)){
  file <- paste0("Facial_recognition/faces/s",n,".pgm")
  pgm_image <- read.pnm(file)
  I_j <- pgm_image@grey
  test <- cbind(test,c(I_j))
}
test <- test[,-1]
X_test <- test - Omega
  
# --- Compute SVD on X^T ---
SVD <- svd(t(X)) 
Theta <- SVD$v  #will have dim pxn as n<<p so R drops extra columns
D <- SVD$d
U <- SVD$u

# --- Plot Scree Plot ---
alpha <- D^2 / (n-1)
plot(seq(1:k), (alpha/cumsum(alpha))[1:k], type='l')

# --- Plot Example Faces ---
plot_face(X[,(56)])
plot_face(X[,(58)])
plot_face(X[,(63)])

# --- Plot Mean Face ---
plot_face(Omega)

# --- Plot Mean Subtracted Example Faces
plot_face(X[,(56)] - Omega)
plot_face(X[,(58)] - Omega)
plot_face(X[,(63)] - Omega)

# --- Plot First 6 Eigen Faces ---
par(mfrow = c(1, 1))
for(h in 1:6){
  plot_face(Theta[,h])
}

# --- Test New Images for All j=1,...,40 ---
id <- c()
MH <- c()
for(j in 1:40){
  #Match test image j to closest training image and display both images
  result <- d_ij(j=j, X=X, X_test=X_test, k=20)
  id[j] <- result$id
  MH[j] <- result$MH
  
  #(Test ID, Matched Training ID, MH dist between them)
  print(paste0("j: ",j, ", Match ", id[j], ", MH: ", MH[j]))
}

# --- Noticeable Cases ---
#10 wrongly matched to 69 MH: 4.56
plot_face(X_test[,10])
plot_face(X[,69])

#36 impressive, matched to 321
plot_face(X_test[,36])
plot_face(X[,321])

# --- Plot MH Values ---
MH_df <- data.frame(Index = 1:40, Value = MH)
ggplot(MH_df, aes(x = Index, y = Value)) +
  geom_point(color = "blue", size = 2) +  
  ggtitle("Scatter Plot of MH[j] Values") + 
  xlab("Index") +  
  ylab("Mahalanobis Distance") +  
  theme_minimal() +
  geom_hline(yintercept = 0, linetype="dashed", color = "grey", size = 0.5) +  
  theme(plot.title = element_text(hjust = 0.5))  +
  geom_point(data = MH_df[10, ], aes(x = Index, y = Value), color = "red", size = 10, shape = 1)



