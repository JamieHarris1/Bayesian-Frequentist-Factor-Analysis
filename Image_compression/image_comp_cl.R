library(jpeg)
library(grid)
library(SpatialPack)
library(abind)

perform_pca <- function(data, q){
  #Image matrices are transposed when read in by jpeg
  X <- t(data)
  
  #Center and scale image
  mu <- colMeans(X)
  sigma <- apply(X, 2, sd)
  X <- (X-mu)/sigma
  
  #Compute SVD
  SVD<- svd(X)
  U <- SVD$u
  D <- diag(SVD$d)
  phi <- SVD$v
  
  #Compute q approx of Z
  Z_unscaled <- U[,1:q] %*% D[1:q,1:q] %*% t(phi[,1:q])
  Z <- (Z_unscaled *sigma) + mu
  
  #Write image to storage
  image <- t(Z)
  return(image)
}



for(q in c(10, 20, 50, 100, 500)){
   #Load in photo
  file_name <- 'cheetah'
  photo <- readJPEG(paste('Image_compression/', file_name,'/CL/cheetah_CL.jpeg',
                          sep=''))
  
  #3 colour channels
  r <- photo[,,1]
  g <- photo[,,2]
  b <- photo[,,3]
  
  photo.r.pca <- perform_pca(r, q)
  photo.g.pca <- perform_pca(g, q)
  photo.b.pca <- perform_pca(b, q)
  image <- list(photo.r.pca, photo.g.pca, photo.b.pca)
  image <- abind(image, along = 3)
  
  writeJPEG(image, paste('Image_compression/', file_name, '/CL/', file_name, '_', q,
                         '_CL_comp.jpeg',sep=''), quality=1)
}
