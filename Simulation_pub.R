####Libraries####
library(data.table)
library(dplyr)
library(fnets)
library(ks)
library(lessR)
library(lubridate)
library(MASS)
library(matrixcalc)
library(mgarch)
library(pracma)
library(reshape2)
library(rmgarch)
library(rockchalk)
library(rugarch)
library(shapes)
library(starvars)
library(StatPerMeCo)
library(tensr)
library(tidyr)
library(vars)
library(xts)


####Simulation####
source('GSI_functions.R')

time_grid <- seq(0, 800, length.out = 1000)  # 1000 time points from 0 to 10
d <- 100     # 20-dimensional covariance matrices
v <- 5     # degrees of freedom = 5 latent processes
lengthscale <- 1.0

set.seed(123)
Sigma <- simulate_wishart_process(time_grid, d, v, lengthscale)
plot.ts(Sigma[1,1,])

####VAR####
##0. VAR on Covariances##
ndays = dim(Sigma)[3]
N = d
window = 700
ntest = 300
eigen1 = matrix(nrow = ndays, ncol = N) 
ntilde = N*(N+1)/2
Sigmavech = matrix(NA, ncol = ntilde, nrow = ndays)
for(l in 1:ndays){
  eigen1[l,] = eigen(Sigma[,,l])$values
  Sigmavech[l,] = matrixcalc::vech(Sigma[,,l])
  print(l)
}
colMeans(eigen1)
pred0 <- matrix(NA, nrow = ntest, ncol = ntilde)
pos0 <- rep(NA, ntest)
for(i in 1:ntest){
  pca_res <- prcomp(Sigmavech[i:(i+window-1),], scale. = TRUE)  # scale. = TRUE if you want to scale the data
  num_pc <- 50
  pc_data <- pca_res$x[, 1:num_pc] 
  var1 <- vars::VAR(pc_data, p = 5, type = 'const')
  foreca0 <- rbindlist(lapply(predict(var1, n.ahead = 1)[[1]], as.data.table))
  loadings <- pca_res$rotation[, 1:num_pc]  # [original dimension x num_pc]
  means <- pca_res$center                   # [1 x original dimension]
  sds <- pca_res$scale                      # [1 x original dimension]
  # Predictions reconstruction
  original_forecast <- foreca0$fcst %*% t(loadings)  # [n.ahead x original dimension]
  original_forecast <- sweep(original_forecast, 2, sds, FUN = "*")
  original_forecast <- sweep(original_forecast, 2, means, FUN = "+")
  pred0[i,] <- original_forecast
  pos0[i] <- isPSD(invvech(pred0[i,]))
  print(i)
}
##1. Approximant##
pred1 <- matrix(NA, nrow = ntest, ncol = ntilde)
pos1 <- rep(NA, ntest)
for(i in 1:ntest){
  pred1[i,] <- vech(approximant(invvech(pred0[i,])))
  pos1[i] <- isPSD(invvech(pred1[i,]))
  print(i)
}

##2.Cholesky##
library(highfrequency)
cholfact <- matrix(NA, nrow = ndays, ncol = ntilde)
tmp <- matrix(nrow = ndays, ncol = ntilde)
for(i in 1:ndays){
  ch1 <- Cholmat(Sigma[,,i],tol = 0.0000000000000000000000000001)
  cholfact[i,] <- t(vech(ch1))
  tmp[i,] <- vech(ch1 %*% t(ch1))
}
pred2 <- matrix(NA, nrow = ntest, ncol = ntilde)
pos2 <- rep(NA, ntest)
for(i in 1:ntest){
  pca_res2 <- prcomp(cholfact[i:(i+window-1),], scale. = TRUE)  # scale. = TRUE if you want to scale the data
  num_pc <- 50
  pc_data2 <- pca_res2$x[, 1:num_pc] 
  var2 <- vars::VAR(pc_data2, p = 5, type = 'const')
  foreca2 <- rbindlist(lapply(predict(var2, n.ahead = 1)[[1]], as.data.table))
  loadings <- pca_res2$rotation[, 1:num_pc]  # [original dimension x num_pc]
  means <- pca_res2$center                   # [1 x original dimension]
  sds <- pca_res2$scale                      # [1 x original dimension]
  # Predictions reconstruction
  original_forecast <- foreca2$fcst %*% t(loadings)  # [n.ahead x original dimension]
  original_forecast <- sweep(original_forecast, 2, sds, FUN = "*")
  original_forecast <- sweep(original_forecast, 2, means, FUN = "+")
  l2 <- t(vech2mat(original_forecast, lowerOnly = T))
  pred2[i,] <- t(vech(crossprod(l2)))
  pos2[i] <- isPSD(invvech(pred2[i,]))
  print(i)
}