library(pixmap)
options(warn = -1)
k <- 200
r <- matrix(NA, nrow=10304)
for(n in 11:409){
  if(!n%%10==0){
    file <- paste0("Facial_recognition/faces/s",n,".pgm")
    pgm_image <- read.pnm(file)
    I_i <- pgm_image@grey
    r <- cbind(r,c(I_i))
  }
}
r <- r[,-1]
Omega <- rowMeans(r)
X <- r - Omega

#test image
test <- matrix(NA, nrow=10304)
for(n in seq(110,4010,100)){
  file <- paste0("Facial_recognition/faces/s",n,".pgm")
  pgm_image <- read.pnm(file)
  I_j <- pgm_image@grey
  test <- cbind(test,c(I_j))
}
test <- test[,-1]
X_test <- test - Omega
  
svd_X <- svd(X) 
Theta <- svd_X$u  #will have dim pxn as n<<p so R drops extra columns
alpha <- svd_X$d

plot(seq(1:k), (alpha/cumsum(alpha))[1:k], type='l')


plot_face <- function(mat){
  plot(pixmapGrey(255* matrix(mat, nrow=112, ncol=92)))
}

d_ij <- function(j,X,X_test, k){
  z_test <- t(Theta[,1:k]) %*% X_test[,j]
  d_i <- c()
  for(i in 1:dim(X)[2]){
    z <- t(Theta[,1:k]) %*% X[,i]
    d_i[i] <- sqrt(t(z-z_test) %*% (z-z_test))
  }
  id <- which.min(d_i)
  plot_face(X[,id])
  plot_face(X_test[,j])
  return(id)
}

#plot mean face
plot_face(Omega)

#plot mean subtracted example faces
plot_face(X[,(56)] - Omega)
plot_face(X[,(58)] - Omega)
plot_face(X[,(63)] - Omega)

#plot eigen faces
par(mfrow = c(1, 1))
for(h in 1:6){
  plot_face(Theta[,h])
}

#reproduce 63,58, 56 and Z values
plot_face(Theta[,1:k] %*% t(Theta[,1:k]) %*% X[,56])
#z values
(t(Theta[,1:k]) %*% X[,56])[1:3]

plot_face(Theta[,1:k] %*% t(Theta[,1:k]) %*% X[,58])
#z values
(t(Theta[,1:k]) %*% X[,58])[1:3]

plot_face(Theta[,1:k] %*% t(Theta[,1:k]) %*% X[,63])
#z values
(t(Theta[,1:k]) %*% X[,63])[1:3]

#test for all j=1,...,40
#prints two images, 1 test and closest image
d_ij(36,X,X_test, 20)
d_ij(10,X,X_test, 20)
#10 wrong
#36 impressive

