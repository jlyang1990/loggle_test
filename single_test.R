##################################################
## installation and preparation
##################################################

## install R package "SINGLE" from pre-packaged binary
# install.packages("SINGLE_1.3.tar.gz", repos=NULL)

## load package "SINGLE"
library(SINGLE)
## load package "loggle"
library(loggle)
## load auxiliary functions for evaluating simulation results, cv model selection, etc.: 
## loggle.eval, loggle.cv.select.full, loggle.sector
source("loggle_test_function.R")


##################################################
## scenario: p=100 time-variant
##################################################

## load data matrix and true precision matrices
load("loggle_data_100.rda")
data <- t(X)

## positions of time points to estimate graphs: 49 equally spaced time points
pos <- round(seq(0.02, 0.98, length=49)*(dim(X)[2]-1)+1)
## lists of h, l1 and l2 for grid search
h.list <- seq(100, 500, 100)  ## bandwidth in kernel smoothed sample covariance matrix
l1.list <- seq(0.1, 0.5, 0.1)  ## tuning parameter of lasso penalty
l2.list <- c(0.5, 0.75, 1, 2, 3)  ## tuning parameter of fused-lasso penalty

## estimate time-varying graphs on grid specified by h.list, l1.list and l2.list
for(h in h.list) {
  
  ## pre-calculate sample covariance matrix with bandwidth h
  ts <- proc.time()
  C <- get_kern_cov(data = data, h = h, kernel = "gaussian")
  te <- proc.time()
  print(sprintf("Time used for get_kern_cov: %.2fs", (te-ts)[3]))
  result_dir <- sprintf("100_h%d/", h)
  save(C, file = sprintf("%sC.RData", result_dir))
  
  ## estimate time-varying graphs on grid specified by l1.list and l2.list
  load(sprintf("%sC.RData", result_dir))
  for(l1 in l1.list) {
    for(l2 in l2.list) {
      ts <- proc.time()
      result <- SINGLE(data = data, C = C, l1 = l1, l2 = l2, verbose = T)
      te <- proc.time()
      print(sprintf("Time used for SINGLE: %.2fs", (te-ts)[3]))
      save(result, file = sprintf("%ssingle_%.2f_%.2f.RData", result_dir, l1, l2))
    }
  }
}

## model selection evaluation for method "single" on grid specified by h.list, l1.list and l2.list
for(h in h.list) {
  for(l1 in l1.list) {
    for(l2 in l2.list) {
      load(sprintf("%ssingle_%.2f_%.2f.RData", result_dir, l1, l2))
      eval.result <- loggle.eval(pos, result$P[pos], Omega.true, X)
      print(sprintf("single method: h %.2f, l1 %.2f, l2 %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f",
                    h, l1, l2, eval.result$FDR[1], eval.result$power[1], eval.result$F1, eval.result$KLdis[1]))
    }
  }
}


##################################################
## scenario: p=500 time-variant
##################################################

## load data matrix and true precision matrices
load("loggle_data_500.rda")
data <- t(X)

## positions of time points to estimate graphs: 49 equally spaced time points
pos <- round(seq(0.02, 0.98, length=49)*(dim(X)[2]-1)+1)
## lists of h, l1 and l2 for grid search
h.list <- c(200, 400)  ## bandwidth in kernel smoothed sample covariance matrix
l1.list <- seq(0.1, 0.5, 0.1)  ## tuning parameter of lasso penalty
l2.list <- c(0.5, 0.75, 1, 2, 3)  ## tuning parameter of fused-lasso penalty

## estimate time-varying graphs on grid specified by h.list, l1.list and l2.list
for(h in h.list) {
  
  ## pre-calculate sample covariance matrix with bandwidth h
  ts <- proc.time()
  C <- get_kern_cov(data = data, h = h, kernel = "gaussian")
  te <- proc.time()
  print(sprintf("Time used for get_kern_cov: %.2fs", (te-ts)[3]))
  result_dir <- sprintf("500_h%d/", h)
  save(C, file = sprintf("%sC.RData", result_dir))
  
  ## estimate time-varying graphs on grid specified by l1.list and l2.list
  load(sprintf("%sC.RData", result_dir))
  for(l1 in l1.list) {
    for(l2 in l2.list) {
      ts <- proc.time()
      result <- SINGLE(data = data, C = C, l1 = l1, l2 = l2, verbose = T)
      te <- proc.time()
      print(sprintf("Time used for SINGLE: %.2fs", (te-ts)[3]))
      save(result, file = sprintf("%ssingle_%.2f_%.2f.RData", result_dir, l1, l2))
    }
  }
}

## model selection evaluation for method "single" on grid specified by h.list, l1.list and l2.list
for(h in h.list) {
  for(l1 in l1.list) {
    for(l2 in l2.list) {
      load(sprintf("%ssingle_%.2f_%.2f.RData", result_dir, l1, l2))
      eval.result <- loggle.eval(pos, result$P[pos], Omega.true, X)
      print(sprintf("single method: h %.2f, l1 %.2f, l2 %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f",
                    h, l1, l2, eval.result$FDR[1], eval.result$power[1], eval.result$F1, eval.result$KLdis[1]))
    }
  }
}


##################################################
## scenario: p=100 time-invariant
##################################################

## load data matrix and true precision matrices
load("loggle_data_100_invar.rda")
data <- t(X)

## positions of time points to estimate graphs: 49 equally spaced time points
pos <- round(seq(0.02, 0.98, length=49)*(dim(X)[2]-1)+1)
## lists of h, l1 and l2 for grid search
h.list <- seq(100, 500, 100)  ## bandwidth in kernel smoothed sample covariance matrix
l1.list <- seq(0.1, 0.5, 0.1)  ## tuning parameter of lasso penalty
l2.list <- c(0.5, 0.75, 1, 2, 3)  ## tuning parameter of fused-lasso penalty

## estimate time-varying graphs on grid specified by h.list, l1.list and l2.list
for(h in h.list) {
  
  ## pre-calculate sample covariance matrix with bandwidth h
  ts <- proc.time()
  C <- get_kern_cov(data = data, h = h, kernel = "gaussian")
  te <- proc.time()
  print(sprintf("Time used for get_kern_cov: %.2fs", (te-ts)[3]))
  result_dir <- sprintf("100_invar_h%d/", h)
  save(C, file = sprintf("%sC.RData", result_dir))
  
  ## estimate time-varying graphs on grid specified by l1.list and l2.list
  load(sprintf("%sC.RData", result_dir))
  for(l1 in l1.list) {
    for(l2 in l2.list) {
      ts <- proc.time()
      result <- SINGLE(data = data, C = C, l1 = l1, l2 = l2, verbose = T)
      te <- proc.time()
      print(sprintf("Time used for SINGLE: %.2fs", (te-ts)[3]))
      save(result, file = sprintf("%ssingle_%.2f_%.2f.RData", result_dir, l1, l2))
    }
  }
}

## model selection evaluation for method "single" on grid specified by h.list, l1.list and l2.list
for(h in h.list) {
  for(l1 in l1.list) {
    for(l2 in l2.list) {
      load(sprintf("%ssingle_%.2f_%.2f.RData", result_dir, l1, l2))
      eval.result <- loggle.eval(pos, result$P[pos], Omega.true, X)
      print(sprintf("single method: h %.2f, l1 %.2f, l2 %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f",
                    h, l1, l2, eval.result$FDR[1], eval.result$power[1], eval.result$F1, eval.result$KLdis[1]))
    }
  }
}