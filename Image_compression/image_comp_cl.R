library(jpeg)
library(abind)
perform_pca <- function(data, k){
  #Image matrices are transposed when read in by jpeg
  X <- t(data)
  
  #Center and scale image
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd)
  X <- (X-mu)
  
  #Compute SVD
  SVD<- svd(X)
  Theta <- SVD$u
  A <- diag(SVD$d)
  U <- SVD$v
  
  #Compute k approx of X
  X_tilde_unscaled <- Theta[,1:k] %*% A[1:k,1:k] %*% t(U[,1:k])
  X_tilde <- (X_tilde_unscaled) + mu
  
  #Write image to storage
  image <- t(X_tilde)
  return(image)
}

for(k in c(10, 20, 50, 100, 500)){
   #Load in photo
  file_name <- 'cheetah'
  photo <- readJPEG(paste('Image_compression/', file_name,'/CL/cheetah_CL.jpeg',
                          sep=''))
  
  #3 colour channels
  r <- photo[,,1]
  g <- photo[,,2]
  b <- photo[,,3]
  
  photo.r.pca <- perform_pca(r, k)
  photo.g.pca <- perform_pca(g, k)
  photo.b.pca <- perform_pca(b, k)
  image <- list(photo.r.pca, photo.g.pca, photo.b.pca)
  image <- abind(image, along = 3)
  
  writeJPEG(image, paste('Image_compression/', file_name, '/CL/', file_name, '_', k,
                         '_CL_comp.jpeg',sep=''), quality=1)
}
