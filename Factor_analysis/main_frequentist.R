library(dplyr)
library(statmod)
library(MASS)
library(psych)
library(lubridate)
library(zoo)
library(readxl)
library(tidyverse)
source("Factor_analysis/Methods/EM.R")
source("Factor_analysis/Methods/PCA.R")
set.seed(123)



# --- DATA ---
# omit the first two rows and the last four rows:
load("Factor_analysis/Macro_data/stockwatson.Rda")
data <- data.frame(stockwatson)
start_date <- as.yearqtr("1959 Q1")
end_date <- as.yearqtr("2008 Q4")
quarterly_dates <- seq(start_date, end_date, by = 1/4)
quarterly_dates <- as.Date(quarterly_dates)
data <- cbind(quarterly_dates, data)
data = data[-c(1:2, 197:200),]
macro_desc <- read_excel("Factor_analysis/Macro_data/macro_desc.xlsx")
colnames(data)[2:dim(data)[2]] <- macro_desc$`Short name`

# --- EDA ---
plot(data$quarterly_dates, data$CPI, type='l',xlab="Date", ylab="Transformed CPI")
longer <- pivot_longer(data[,c(1,4,9,10)], cols =2:4 , names_to = "group", values_to = "values")

ggplot(longer, aes(x = quarterly_dates, y = values, color = group)) +
  geom_line() +  # Add lines for each time series
  labs(title = "Multiple Time Series Plot", x = "Time", y = "Value") +  # Set plot labels
  theme_minimal() +  # Set a minimal theme
  scale_color_manual(values = c("blue", "orange", "darkgreen"))  # Customize color scale


k <- 25
# --- PCA FA ---
n <- dim(data[,-1])[1]
p <- dim(data[,-1])[2]
pca_result <- PCA(data[,-1])
theta <- pca_result$theta[,1:k]
D2 <- pca_result$D2[1:k]
lambda <- matvec(theta, as.numeric(lapply(D2, sqrt)))
lambda <- as.numeric(varimax(lambda)$loadings)
lambda <- matrix(lambda, nrow=p, ncol=k)
rownames(lambda) <- colnames(data[,-1])

sum(D2[1:2]) / sum(D2)
norm(cov(data[,-1]) - lambda%*%t(lambda), type = "F")

# --- EM FA ---
result <- em(data[,-1], k=k, iter=200)

lambda <- result$lambda
sigma <- result$sigma

lambda <- as.numeric(varimax(lambda)$loadings)
lambda <- matrix(lambda, nrow=p, ncol=k)
rownames(lambda) <- colnames(data[,-1])
View(lambda)

norm(cov(data[,-1]) - lambda%*%t(lambda), type = "F")

