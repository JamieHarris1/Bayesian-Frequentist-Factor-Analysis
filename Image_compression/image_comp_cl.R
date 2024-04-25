library(jpeg)
library(abind)

perform_pca <- function(X, k){
  X_bar <- colMeans(X)
  X <- (X-X_bar)
  
  # --- Compute SVD ---
  SVD <- svd(X)
  Theta <- SVD$v
  D <- diag(SVD$d)
  U <- SVD$u
  
  # --- Low Rank SVD Approx ---
  X_k_SVD <- Theta[,1:k] %*% D[1:k,1:k] %*% t(U[,1:k])
  image <- t(X_k_SVD)  + X_bar
  return(image)
}

# --- Load in photo ---
file_name <- 'cheetah'
X <- readJPEG(paste('Image_compression/', file_name,'/CL/cheetah_CL.jpeg',
                        sep=''))

# --- Split 3 Colour Channels ---
r <- X[,,1]
g <- X[,,2]
b <- X[,,3]

# --- Image Compression ---
for(k in c(10, 20, 50, 100, 500)){
  photo.r.pca <- perform_pca(r, k)
  photo.g.pca <- perform_pca(g, k)
  photo.b.pca <- perform_pca(b, k)
  image <- list(photo.r.pca, photo.g.pca, photo.b.pca)
  image <- abind(image, along = 3)
  writeJPEG(image, paste('Image_compression/', file_name, '/CL/', file_name, '_', k,
                         '_CL_comp.jpeg',sep=''), quality=1)
}
