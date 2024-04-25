# -- Gamma pdf plot for MH within step --
shape <- 10
scale <- 0.5

x <- seq(0, 10, length.out = 100)
pdf_values <- dgamma(x, shape=shape, scale=scale)
plot(x, pdf_values, type = "l", 
     xlab = "x", ylab = "f(x)", 
     col = "black", lwd = 2)





