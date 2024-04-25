library(jpeg)

# --- Load in photo ---
file_name <- 'cheetah'
X <- readJPEG(paste('Image_compression/', file_name,'/BW/cheetah_BW.jpeg',
                        sep=''))
X_bar <- colMeans(X)
X <- (X-X_bar)

# --- Compute SVD ---
SVD <- svd(X)
Theta <- SVD$v
D <- diag(SVD$d)
U <- SVD$u

# --- Image Compression ---
for(k in c(10, 20, 50, 100, 500)){
  X_k_SVD <- Theta[,1:k] %*% D[1:k,1:k] %*% t(U[,1:k])
  image <- t(X_k_SVD) + X_bar
  writeJPEG(image, paste('Image_compression/', file_name, '/BW/', file_name,"_",
                         k, '_BW_comp.jpeg',sep=''), quality=1)
  
}

