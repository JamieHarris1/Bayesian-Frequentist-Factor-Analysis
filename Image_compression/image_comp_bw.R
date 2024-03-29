library(jpeg)
library(grid)
library(SpatialPack)



#Load in photo
file_name <- 'cheetah'
BW <- readJPEG(paste('Image_compression/', file_name,'/BW/cheetah_BW.jpeg',
                        sep=''))

#Image matrices are transposed when read in by jpeg
X <- t(BW)

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
q=c(10, 20, 50, 100, 500)
for(i in q){
  Z_unscaled <- U[,1:i] %*% D[1:i,1:i] %*% t(phi[,1:i])
  Z <- (Z_unscaled *sigma) + mu
  
  
  #Write image to storage
  image <- t(Z)
  writeJPEG(image, paste('Image_compression/', file_name, '/BW/', file_name,"_",
                         i, '_BW_comp.jpeg',sep=''), quality=1)
  
}

