library(jpeg)

#Load in photo
file_name <- 'cheetah'
BW <- readJPEG(paste('Image_compression/', file_name,'/BW/cheetah_BW.jpeg',
                        sep=''))

#Image matrices are transposed when read in by jpeg
X <- t(BW)

#Center and scale image
mu <- colMeans(X)
X <- (X-mu)

#Compute SVD
SVD<- svd(X)
Theta <- SVD$u
A <- diag(SVD$d)
U <- SVD$v


#Compute k approx of X
for(k in c(10, 20, 50, 100, 500)){
  X_tilde_unscaled <- Theta[,1:k] %*% A[1:k,1:k] %*% t(U[,1:k])
  X_tilde <- (X_tilde_unscaled) + mu
  
  
  #Write image to storage
  image <- t(X_tilde)
  writeJPEG(image, paste('Image_compression/', file_name, '/BW/', file_name,"_",
                         k, '_BW_comp.jpeg',sep=''), quality=1)
  
}

