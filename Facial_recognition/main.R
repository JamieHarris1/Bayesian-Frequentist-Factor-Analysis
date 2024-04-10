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
v <- svd_X$v

plot(seq(1:k), (alpha/cumsum(alpha))[1:k], type='l')


plot_face <- function(mat){
  plot(pixmapGrey(255* matrix(mat, nrow=112, ncol=92)))
}

d_ij <- function(j,X,X_test, k){
  z_test <- t(Theta[,1:k]) %*% X_test[,j]
  d_i <- c()
  for(i in 1:dim(X)[2]){
    z <- t(Theta[,1:k]) %*% X[,i]
    S <- cov(z,z_test)
    d_i[i] <- sqrt(1/S * t(z-z_test) %*% (z-z_test))
  }
  id <- which.min(d_i)
  MH <- d_i[id]
  plot_face(X[,id])
  plot_face(X_test[,j])
  return(list(id=id, MH=MH))
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

#test for all j=1,...,40
#prints two images, 1 test and closest image

id <- c()
MH <- c()
for(j in 1:40){
  result <- d_ij(j=j,X=X,X_test=X_test, k=20)
  id[j] <- result$id
  MH[j] <- result$MH
  print(paste0("j: ",j, ", Match ", id[j], ", MH: ", MH[j]))
}


#10 wrongly matched to 69 MH: 4.56
plot_face(X_test[,10])
plot_face(X[,69])

#36 impressive, matched to 321
plot_face(X_test[,36])
plot_face(X[,321])

#plot MH values
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



