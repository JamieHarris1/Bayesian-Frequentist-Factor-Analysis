library(dplyr)
library(statmod)
library(MASS)
library(psych)
library(lubridate)
library(zoo)
library(readxl)
library(tidyverse)
library(pheatmap)
source("Factor_analysis/Methods/EM.R")
source("Factor_analysis/Methods/PCA.R")
set.seed(0)

# --- DATA ---
load("Factor_analysis/Macro_data/stockwatson.Rda")
data <- data.frame(stockwatson)
start_date <- as.yearqtr("1959 Q1")
end_date <- as.yearqtr("2008 Q4")
quarterly_dates <- seq(start_date, end_date, by = 1/4)
quarterly_dates <- as.Date(quarterly_dates)
data <- cbind(quarterly_dates, data)
# omit the first two rows and the last four rows:
data = data[-c(1:2, 197:200),]
macro_desc <- read_excel("Factor_analysis/Macro_data/macro_desc.xlsx")
colnames(data)[2:dim(data)[2]] <- macro_desc$`Short name`
data[,-1] <- scale(data[,-1], center = TRUE)
X <- t(data[,-1])

# --- EDA ---
plot(data$quarterly_dates, data$CPI, type='l',xlab="Date", ylab="Transformed CPI")
longer <- pivot_longer(data[,c(1,4,9,10)], cols =2:4 , names_to = "group", values_to = "values")

ggplot(longer, aes(x = quarterly_dates, y = values, color = group)) +
  geom_line() +  # Add lines for each time series
  labs(title = "Multiple Time Series Plot", x = "Time", y = "Value") +  # Set plot labels
  theme_minimal() +  # Set a minimal theme
  scale_color_manual(values = c("blue", "orange", "darkgreen"))  # Customize color scale

pheatmap(cor(data[,-1]), 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Covariance Matrix Heatmap",
         display_numbers = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)



# --- PCA FA ---
k <- 25
p <- dim(X)[1]
n <- dim(X)[2]

pca_result <- PCA(X)
Theta_k <- pca_result$Theta[,1:k]
alpha <- pca_result$alpha
alpha_k <- alpha[1:k]
Lambda <- matvec(Theta_k, as.numeric(lapply(alpha_k, sqrt)))
Lambda <- as.numeric(varimax(Lambda)$loadings)
Lambda <- matrix(Lambda, nrow=p, ncol=k)
rownames(Lambda) <- colnames(data[,-1])
View(Lambda)

pheatmap(Lambda, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Covariance Matrix Heatmap",
         display_numbers = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

#variance explained by k=6 and k=25
print(paste0("k=6: ", sum(alpha[1:6]) / sum(alpha)))
print(paste0("k=25: ", sum(alpha[1:25]) / sum(alpha)))

#Frobenious of error 
norm(cov(t(X)) - Lambda%*%t(Lambda), type = "F")

# --- EM FA ---
result <- em(X, k=k, iter=50000)

Lambda <- result$Lambda
Sigma <- result$Sigma

Lambda <- as.numeric(varimax(Lambda)$loadings)
Lambda <- matrix(Lambda, nrow=p, ncol=k)
rownames(Lambda) <- colnames(data[,-1])
View(Lambda)

#Frobenious of error 
norm(cov(t(X)) - Lambda%*%t(Lambda), type = "F")

pheatmap(Lambda, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Covariance Matrix Heatmap",
         display_numbers = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
