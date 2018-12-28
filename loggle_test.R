##################################################
## installation and preparation
##################################################

## install package dependencies (R with version later than 3.0.2 needed)
# install.packages(c("Matrix", "doParallel", "igraph", "glasso", "sm"))

## install R package "loggle" from pre-packaged binary
## in Linux/Mac OS X:
# install.packages("loggle_1.0.tar.gz", repos=NULL)
## in Windows:
# install.packages("loggle_1.0.zip", repos=NULL)

## alternative: install R package "loggle" from source files in GitHub repository (R package "devtools" needed)
# install.packages("devtools")
# library(devtools)
# install_github(repo="jlyang1990/loggle")

## load package "loggle"
library(loggle)
## help of the main function "loggle"
?loggle


##################################################
## an example
##################################################
data(example)  ## load example dataset
## data matrix and true precision matrices
X <- example$X
Omega.true <- example$Omega.true
dim(X)  ## dimension of data matrix
p <- nrow(X)  ## number of variables

## positions of time points to estimate graphs
pos <- round(seq(0.02, 0.98, length=25)*(ncol(X)-1)+1)
K <- length(pos)
## estimate time-varying graphs and conduct model selection via cross-validation
## num.thread can be set as large as number of cores on a multi-core machine 
## (however when p is large, memory overflow should also be taken caution of)
ts <- proc.time()
result <- loggle.cv(X, pos, h.list = c(0.15, 0.2, 0.25, 0.3), d.list = c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 1), 
                    lambda.list = c(0.15, 0.2, 0.25, 0.3), fit.type = "pseudo", cv.vote.thres = 1, num.thread = 1)
te <- proc.time()
sprintf("Time used for loggle.cv: %.2fs", (te-ts)[3])

## optimal values of h, d and lambda, and number of selected edges at each time point
select.result <- result$cv.select.result
print(cbind("time" = seq(0.02, 0.98, length=25), "h.opt" = rep(select.result$h.opt, K), "d.opt" = select.result$d.opt,
            "lambda.opt" = select.result$lambda.opt, "edge.num.opt" = select.result$edge.num.opt))

## false discovery rate (FDR) and power based on true precision matrices for selected model
edge.num.opt <- select.result$edge.num.opt
edge.num.true <- sapply(1:K, function(i) (sum(Omega.true[[pos[i]]]!=0)-p)/2)
edge.num.overlap <- sapply(1:K, function(i) (sum(select.result$adj.mat.opt[[i]] & Omega.true[[pos[i]]])-p)/2)
perform.matrix <- cbind("FDR" = 1 - edge.num.overlap / edge.num.opt, "power" = edge.num.overlap / edge.num.true)
print(apply(perform.matrix, 2, mean))


##################################################
## simulation in the paper:
## note: we recommend using a server to conduct the following simulations
##################################################
## load auxiliary functions for evaluating simulation results, cv model selection, etc.: 
## loggle.eval, loggle.cv.select.full, loggle.sector
source("loggle_test_function.R")


##################################################
## scenario: p=100 time-variant
##################################################

## load data matrix, true precision matrices and pre-tuned hyper-parameters(h, d and lambda)
load("loggle_data_100.rda")

## positions of time points to estimate graphs: 49 equally spaced time points
pos <- round(seq(0.02, 0.98, length=49)*(dim(X)[2]-1)+1)
## lists of h, d and lambda for grid search
h.list <- seq(0.1, 0.3, 0.05)  ## bandwidth in kernel smoothed sample covariance/correlation matrix
d.list <- c(0, 0.001, 0.01, seq(0.025, 0.1, 0.025), seq(0.15, 0.3, 0.05), 1)  ## width of neighborhood
lambda.list <- seq(0.15, 0.35, length = 11)  ## tuning parameter of lasso penalty

## fit on a server (72 cores and 256 RAM): default number of thread = 25
## time cost: 3750s - recommended using server rather than laptop

## estimate time-varying graphs and conduct model selection via cross-validation
ts <- proc.time()
result <- loggle.cv(X, pos, h.list, d.list, lambda.list, early.stop.thres = 3, fit.type = "pseudo", 
                    epi.abs = c(rep(1e-5, 11), 1e-4), epi.rel = c(rep(1e-3, 11), 1e-2), detrend = FALSE, 
                    num.thread = 25)
te <- proc.time()
sprintf("Time used for loggle.cv: %.2fs", (te-ts)[3])
  
## model selection evaluation for method "loggle"
select.result <- loggle.cv.select.full(result, select.type = "all_flexible", cv.vote.thres = 1, method = "loggle")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("loggle method: h.opt %.2f, d.opt %.2f(%.2f), lambda.opt %.2f(%.2f), FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f",
        select.result$h.opt, mean(select.result$d.opt), sd(select.result$d.opt), mean(select.result$lambda.opt), 
        sd(select.result$lambda.opt), eval.result$FDR[1], eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# loggle method: h.opt 0.15, d.opt 0.27(0.03), lambda.opt 0.19(0.02), FDR 0.196, power 0.702, F1 0.747, KLdis 2.284
  
## model selection evaluation for method "kernel"
select.result <- loggle.cv.select.full(result, select.type = "all_fixed", cv.vote.thres = 1, method = "kernel")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("kernel method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f", 
        select.result$h.opt, mean(select.result$d.opt), mean(select.result$lambda.opt), eval.result$FDR[1], 
        eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# kernel method: h.opt 0.20, d.opt 0.00, lambda.opt 0.21, FDR 0.063, power 0.571, F1 0.703, KLdis 2.690
  
## model selection evaluation for method "invariant"
select.result <- loggle.cv.select.full(result, select.type = "all_fixed", cv.vote.thres = 1, method = "invariant")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("invariant method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f", 
        select.result$h.opt, mean(select.result$d.opt), mean(select.result$lambda.opt), eval.result$FDR[1], 
        eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# invariant method: h.opt 0.15, d.opt 1.00, lambda.opt 0.15, FDR 0.583, power 0.678, F1 0.514, KLdis 2.565

## time-varying graphs from model selection can be recovered using pre-tuned hyperparameters(h, d and lambda)
## this can be fitted on a laptop: default number of thread = 1

## model selection evaluation for method "loggle" using pre-tuned hyperparameters
## method "loggle": possibly different tuning parameters for different positions
## time cost: num.thread=4: 89.10s, num.thread=2: 158.80s, num.thread=1: 296.00s
ts <- proc.time()
result <- loggle.cv.vote(X, pos, h = tuned.parameter$loggle$h, d = tuned.parameter$loggle$d, 
                         lambda = tuned.parameter$loggle$lambda, fit.type = "pseudo", cv.vote.thres = 1, 
                         epi.abs = 1e-5, epi.rel = 1e-3, detrend = FALSE, num.thread = 1)
te <- proc.time()
sprintf("Time used for loggle: %.2fs", (te-ts)[3])
eval.result <- loggle.eval(pos, result$Omega, Omega.true, X)
sprintf("loggle method: h.opt %.2f, d.opt %.2f(%.2f), lambda.opt %.2f(%.2f), FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f",
        tuned.parameter$loggle$h, mean(tuned.parameter$loggle$d), sd(tuned.parameter$loggle$d), 
        mean(tuned.parameter$loggle$lambda), sd(tuned.parameter$loggle$lambda), eval.result$FDR[1], 
        eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# loggle method: h.opt 0.15, d.opt 0.27(0.03), lambda.opt 0.19(0.02), FDR 0.196, power 0.702, F1 0.747, KLdis 2.285

## model selection evaluation for method "kernel" using pre-tuned hyperparameters
## method "kernel": fitting one model for each position independently without structure smoothing
## time cost: num.thread=4: 19.72s, num.thread=2: 21.00s, num.thread=1: 24.70s
ts <- proc.time()
result <- loggle.cv.vote(X, pos, h = tuned.parameter$kernel$h, d = tuned.parameter$kernel$d, 
                         lambda = tuned.parameter$kernel$lambda, fit.type = "pseudo", cv.vote.thres = 1, 
                         epi.abs = 1e-5, epi.rel = 1e-3, detrend = FALSE, num.thread = 1)
te <- proc.time()
sprintf("Time used for kernel: %.2fs", (te-ts)[3])
eval.result <- loggle.eval(pos, result$Omega, Omega.true, X)
sprintf("kernel method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f", 
        tuned.parameter$kernel$h, mean(tuned.parameter$kernel$d), mean(tuned.parameter$kernel$lambda), 
        eval.result$FDR[1], eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# kernel method: h.opt 0.20, d.opt 0.00, lambda.opt 0.21, FDR 0.063, power 0.571, F1 0.703, KLdis 2.691

## model selection evaluation for method "invariant" using pre-tuned hyperparameters
## method "invariant": fitting time-invariant graphs with same structures for all positions
## time cost: num.thread=1: 37.33s
ts <- proc.time()
result <- loggle.cv.vote(X, pos, h = tuned.parameter$invar$h, d = tuned.parameter$invar$d, 
                         lambda = tuned.parameter$invar$lambda, fit.type = "pseudo", cv.vote.thres = 1, 
                         epi.abs = 1e-4, epi.rel = 1e-2, detrend = FALSE, num.thread = 1)
te <- proc.time()
sprintf("Time used for invariant: %.2fs", (te-ts)[3])
eval.result <- loggle.eval(pos, result$Omega, Omega.true, X)
sprintf("invariant method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f", 
        tuned.parameter$invar$h, mean(tuned.parameter$invar$d), mean(tuned.parameter$invar$lambda), 
        eval.result$FDR[1], eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# invariant method: h.opt 0.15, d.opt 1.00, lambda.opt 0.15, FDR 0.584, power 0.678, F1 0.514, KLdis 2.570


##################################################
## scenario: p=500 time-variant
##################################################

## load data matrix and true precision matrices
load("loggle_data_500.rda")

## positions of time points to estimate graphs: 49 equally spaced time points
pos <- round(seq(0.02, 0.98, length=49)*(dim(X)[2]-1)+1)
## lists of h, d and lambda for grid search
h.list <- c(0.15, 0.2)
d.list <- c(0, 0.001, 0.01, seq(0.025, 0.1, 0.025), seq(0.15, 0.3, 0.05), 1)
lambda.list <- seq(0.15, 0.35, length = 11)

## fit on a server (72 cores and 256 RAM): default number of thread = 25
## time cost: 85390s - recommended using server rather than laptop
## note: requiring too many cores may cause memory overflow and the following error message 
##########
# Error in Omega.list[, -D, k] <- result.k$Omega.list :
# number of items to replace is not a multiple of replacement length
##########

## estimate time-varying graphs and conduct model selection via cross-validation
ts <- proc.time()
result <- loggle.cv(X, pos, h.list, d.list, lambda.list, early.stop.thres = 1, fit.type = "pseudo", 
                    epi.abs = 1e-4, epi.rel = 1e-2, detrend = FALSE, num.thread = 25)
te <- proc.time()
sprintf("Time used for loggle.cv: %.2fs", (te-ts)[3])
  
## model selection evaluation for method "loggle"
select.result <- loggle.cv.select.full(result, select.type = "all_flexible", method = "loggle")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("loggle method: h.opt %.2f, d.opt %.2f(%.2f), lambda.opt %.2f(%.2f), FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f",
        select.result$h.opt, mean(select.result$d.opt), sd(select.result$d.opt), mean(select.result$lambda.opt), 
        sd(select.result$lambda.opt), eval.result$FDR[1], eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# loggle method: h.opt 0.15, d.opt 0.29(0.02), lambda.opt 0.21(0.02), FDR 0.215, power 0.613, F1 0.678, KLdis 9.564
  
## model selection evaluation for method "kernel"
select.result <- loggle.cv.select.full(result, select.type = "all_fixed", method = "kernel")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("kernel method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f", 
        select.result$h.opt, mean(select.result$d.opt), mean(select.result$lambda.opt), eval.result$FDR[1], 
        eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# kernel method: h.opt 0.20, d.opt 0.00, lambda.opt 0.27, FDR 0.035, power 0.399, F1 0.561, KLdis 11.818
  
## model selection evaluation for method "invariant"
select.result <- loggle.cv.select.full(result, select.type = "all_fixed", method = "invariant")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("invariant method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f", 
        select.result$h.opt, mean(select.result$d.opt), mean(select.result$lambda.opt), eval.result$FDR[1], 
        eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# invariant method: h.opt 0.15, d.opt 1.00, lambda.opt 0.17, FDR 0.590, power 0.597, F1 0.478, KLdis 10.608


##################################################
## scenario: p=100 time-invariant
##################################################

## load data matrix and true precision matrices
load("loggle_data_100_invar.rda")

## positions of time points to estimate graphs: 49 equally spaced time points
pos <- round(seq(0.02, 0.98, length=49)*(dim(X)[2]-1)+1)
## lists of h, d and lambda for grid search
h.list <- seq(0.1, 0.3, 0.05)
d.list <- c(0, 0.001, 0.01, seq(0.025, 0.1, 0.025), seq(0.15, 0.3, 0.05), 1)
lambda.list <- seq(0.15, 0.35, length = 11)

## fit on a server (72 cores and 256 RAM): default number of thread = 25
## time cost: 9669s - recommended using server rather than laptop

## estimate time-varying graphs and conduct model selection via cross-validation
ts <- proc.time()
result <- loggle.cv(X, pos, h.list, d.list, lambda.list, early.stop.thres = 1.25, fit.type = "pseudo", 
                    epi.abs = c(rep(1e-5, 11), 1e-4), epi.rel = c(rep(1e-3, 11), 1e-2), detrend = FALSE, 
                    num.thread = 25)
te <- proc.time()
sprintf("Time used for loggle.cv: %.2fs", (te-ts)[3])
  
## model selection evaluation for method "loggle"
select.result <- loggle.cv.select.full(result, select.type = "all_flexible", cv.vote.thres = 1, method = "loggle")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("loggle method: h.opt %.2f, d.opt %.2f(%.2f), lambda.opt %.2f(%.2f), FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f",
        select.result$h.opt, mean(select.result$d.opt), sd(select.result$d.opt), mean(select.result$lambda.opt), 
        sd(select.result$lambda.opt), eval.result$FDR[1], eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# loggle method: h.opt 0.10, d.opt 0.96(0.17), lambda.opt 0.20(0.02), FDR 0.000, power 0.978, F1 0.988, KLdis 1.559
  
## model selection evaluation for method "kernel"
select.result <- loggle.cv.select.full(result, select.type = "all_fixed", cv.vote.thres = 1, method = "kernel")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("kernel method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f", 
        select.result$h.opt, mean(select.result$d.opt), mean(select.result$lambda.opt), eval.result$FDR[1], 
        eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# kernel method: h.opt 0.15, d.opt 0.00, lambda.opt 0.23, FDR 0.042, power 0.509, F1 0.598, KLdis 3.168
  
## model selection evaluation for method "invariant"
select.result <- loggle.cv.select.full(result, select.type = "all_fixed", cv.vote.thres = 1, method = "invariant")
eval.result <- loggle.eval(pos, select.result$adj.mat.opt, Omega.true, X)
sprintf("invariant method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, FDR %.3f, power %.3f, F1 %.3f, KLdis %.3f", 
        select.result$h.opt, mean(select.result$d.opt), mean(select.result$lambda.opt), eval.result$FDR[1], 
        eval.result$power[1], eval.result$F1, eval.result$KLdis[1])
# invariant method: h.opt 0.10, d.opt 1.00, lambda.opt 0.19, FDR 0.000, power 1.000, F1 1.000, KLdis 1.531


##################################################
## real data application
##################################################

## load stock price dataset
data(stockdata)
date.index <- rownames(stockdata$data)
stock.sector <- stockdata$info[, "sector"]

## select stocks from sectors Information Technology, Consumer Discretionary, Consumer Staples, Financials, Industrials
## 283 stocks, 1008 time points (trading days)
sp.it <- t(stockdata$data[date.index < "2011-01-01", stock.sector == "Information Technology"])
sp.cd <- t(stockdata$data[date.index < "2011-01-01", stock.sector == "Consumer Discretionary"])
sp.cs <- t(stockdata$data[date.index < "2011-01-01", stock.sector == "Consumer Staples"])
sp.f <- t(stockdata$data[date.index < "2011-01-01", stock.sector == "Financials"])
sp.i <- t(stockdata$data[date.index < "2011-01-01", stock.sector == "Industrials"])
sp <- rbind(sp.it, sp.cd, sp.cs, sp.f, sp.i)
sp.num <- c(nrow(sp.it), nrow(sp.cd), nrow(sp.cs), nrow(sp.f), nrow(sp.i))
sector.num <- length(sp.num)

## construct data matrix by taking log ratio of prices between adjacent time points
p <- dim(sp)[1]
N <- dim(sp)[2]-1
X <- matrix(0, p, N)
for(i in 1:p) {
  X[i, ] <- scale(log(sp[i, -1] / sp[i, -(N+1)]))
}
dim(X)  ## dimension of data matrix
# 283 1007

## positions of time points to estimate graphs: 201 equally spaced time points
K <- 201
pos <- round(seq(0.005, 0.995, length=K)*(ncol(X)-1)+1)
loc <- seq(0.005, 0.995, length=K)
## lists of h, d and lambda for grid search
h.list <- c(0.1, 0.15)
d.list <- c(0, 0.001, 0.01, seq(0.025, 0.1, 0.025), seq(0.15, 0.3, 0.05), 1)
lambda.list <- 10 ^ seq(-2, 0, length = 21)

## fit on a server (72 cores and 256 RAM): default number of thread = 23
## recommend using server rather than laptop
## warning: time consuming (a few days per "h" for fitting 201 time points) 
## you can skip the next block and directly go to "model selection results" where pre-selected models are stored 
## such that you do not need conduct grid search and can directly examine the results for selected models instead

##########
## estimate time-varying graphs and conduct model selection via cross-validation
ts <- proc.time()
result <- loggle.cv(X, pos, h.list, d.list, lambda.list, early.stop.thres = 8, fit.type = "pseudo", 
                    epi.abs = 1e-4, epi.rel = 1e-2, num.thread = 23)
te <- proc.time()
sprintf("Time used for loggle.cv: %.2fs", (te-ts)[3])
## estimate time-varying graphs and conduct model selection via cross-validation for method "kernel"
ts <- proc.time()
result.kernel <- loggle.cv(X, pos, h.list, d.list = 0, lambda.list, early.stop.thres = 15, fit.type = "pseudo", 
                            epi.abs = 1e-4, epi.rel = 1e-2, num.thread = 23)
te <- proc.time()
sprintf("Time used for loggle.cv for method kernel: %.2fs", (te-ts)[3])
  
## model selection result for method "loggle"
select.result.loggle <- loggle.cv.select.full(result, select.type = "all_flexible", method = "loggle")
sprintf("loggle method: h.opt %.2f, d.opt %.2f(%.2f), lambda.opt %.2f(%.2f), edge.num.opt %.1f(%.1f), cv.score %.2f", 
        select.result.loggle$h.opt, mean(select.result.loggle$d.opt), sd(select.result.loggle$d.opt), 
        mean(select.result.loggle$lambda.opt), sd(select.result.loggle$lambda.opt), 
        mean(select.result.loggle$edge.num.opt), sd(select.result.loggle$edge.num.opt), 
        select.result.loggle$cv.score.opt)
# loggle method: h.opt 0.10, d.opt 0.56(0.43), lambda.opt 0.52(0.16), edge.num.opt 819.4(331.0), cv.score 123.06
  
## model selection result for method "kernel"
select.result.kernel <- loggle.cv.select.full(result.kernel, select.type = "all_fixed", method = "kernel")
sprintf("kernel method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, edge.num.opt %.1f(%.1f), cv.score %.2f", 
        select.result.kernel$h.opt, mean(select.result.kernel$d.opt), mean(select.result.kernel$lambda.opt), 
        mean(select.result.loggle$edge.num.opt), sd(select.result.loggle$edge.num.opt), 
        select.result.kernel$cv.score.opt)
# kernel method: h.opt 0.15, d.opt 0.00, lambda.opt 0.16, edge.num.opt 1103.5(487.1), cv.score 160.14
  
## model selection result for method "invariant"
select.result.invar <- loggle.cv.select.full(result, select.type = "all_fixed", method = "invariant")
sprintf("invariant method: h.opt %.2f, d.opt %.2f, lambda.opt %.2f, edge.num.opt %.1f(%.1f), cv.score %.2f", 
        select.result.invar$h.opt, mean(select.result.invar$d.opt), mean(select.result.invar$lambda.opt), 
        mean(select.result.loggle$edge.num.opt), sd(select.result.loggle$edge.num.opt), 
        select.result.invar$cv.score.opt)
# invariant method: h.opt 0.10, d.opt 1.00, lambda.opt 0.50, edge.num.opt 811.0(0.0), cv.score 130.68
##########

## model selection results
## load pre-stored results from loggle.cv.select.full
load("loggle_select_result_real.rda")

## construct summary matrix containing number of within-sector and cross-sector edges
edge.sector.loggle <- loggle.sector(select.result.loggle, sp.num, h = 0.04)
edge.sector.kernel <- loggle.sector(select.result.kernel, sp.num, h = 0.04)
edge.sector.invar <- loggle.sector(select.result.invar, sp.num, h = 0.04)

## number of edges for methods loggle, kernel and invariant
sector.plot.loggle <- edge.sector.loggle$edge.sector.sm
sector.plot.kernel <- edge.sector.kernel$edge.sector.sm
sector.plot.invar <- edge.sector.invar$edge.sector.sm
plot(loc, sector.plot.loggle[sector.num+2, ], type = "l", lwd = 2, ylim = c(0, 2000), xlab = "", ylab = "# of edges", 
     xaxt = "n", cex.lab = 1.2)
lines(loc, sector.plot.kernel[sector.num+2, ], type = "l", lty = 2, lwd = 2)
lines(loc, sector.plot.invar[sector.num+2, ], type = "l", lty = 3, lwd = 2)
legend("topright", c("loggle", "kernel", "invar"), lty = 1:3, lwd = 2, cex = 1)
axis(side = 1, seq(0.005, 0.995, length=11), labels = date.index[round(seq(0.005, 0.995, length=11)*(N-1)+1)], las = 2, 
     cex.axis = 0.8)

## proportion of within-sector edges (among total number of detected edges) for methods loggle, kernel and invariant 
plot(loc, colSums(sector.plot.loggle[1:sector.num, ])/sector.plot.loggle[sector.num+2, ], type = "l", lwd = 2, 
     ylim = c(0.2, 1), xlab = "", ylab = "proportion of within-sector edges", xaxt = "n", cex.lab = 1.2)
lines(loc, colSums(sector.plot.kernel[1:sector.num, ])/sector.plot.kernel[sector.num+2, ], type = "l", lty = 2, lwd = 2)
lines(loc, colSums(sector.plot.invar[1:sector.num, ])/sector.plot.invar[sector.num+2, ], type = "l", lty = 3, lwd = 2)
legend("topright", c("loggle", "kernel", "invar"), lty = 1:3, lwd = 2, cex = 1)
axis(side = 1, seq(0.005, 0.995, length=11), labels = date.index[round(seq(0.005, 0.995, length=11)*(N-1)+1)], las = 2, 
     cex.axis = 0.8)

## sector-wise percentage of presence of within-sector edges (defined as the ratio between the number of detected 
## within-sector edges and the total number of possible within-sector edges for a given sector)
## and percentage of presence of cross-sector edges (defined as the ratio between the number of detected cross-sector
## edges and the total number of possible cross-sector edges) of "loggle" fitted graphs
sector.plot <- edge.sector.loggle$edge.sector.sm
plot(loc, sector.plot[sector.num+1, ]/(p*(p-1)/2-sum(sp.num*(sp.num-1)/2)), type = "l", lwd = 2, col = 1, 
     ylim = c(0, 0.2), xlab = "", ylab = "percentage of edge presence", xaxt = "n", cex.lab = 1.2)
for(i in 1:sector.num){
  lines(loc, sector.plot[i, ]/(sp.num[i]*(sp.num[i]-1)/2), col = rainbow(sector.num)[i], lwd = 2)
}
legend("topright", c("I. T.", "Cons. Disc.", "Cons. Staples", "Financials", "Industrials", "Cross Sectors"), lty = 1, 
       lwd = 2, col = c(rainbow(sector.num), 1), cex = 0.8)
axis(side = 1, seq(0.005, 0.995, length=11), labels = date.index[round(seq(0.005, 0.995, length=11)*(N-1)+1)], las = 2,
     cex.axis = 0.8)

## "loggle" fitted graphs at five time points
library(igraph)
pos.plot <- pos[round(seq(5, length(pos)-4, length = 5))]
for(k in 1:length(pos.plot)) {
  adj.matrix <- select.result.loggle$adj.mat.opt[[which(pos == pos.plot[k])]]
  net <- graph.adjacency(adj.matrix, mode = "undirected", diag = FALSE)
  set.seed(0)
  V(net)$color <- unlist(sapply(1:sector.num, function(i) rep(rainbow(sector.num)[i], sp.num[i])))
  E(net)$color <- "darkgray"
  E(net)$width <- 0.6
  plot(net, vertex.size = 3, vertex.label = NA, layout = layout.fruchterman.reingold)
  title(date.index[pos.plot[k]], cex.main = 2)
  if(k == length(pos.plot)) {
    legend(0.45, -0.35, c("I. T.", "Cons. Disc.", "Cons. Staples", "Financials", "Industrials"), pch = 21, 
           col = "black", pt.bg = rainbow(sector.num), cex = 1.2, y.intersp = 0.7, bty = 'n')
  }
}

## time-varying graphs from model selection can be recovered using pre-tuned hyperparameters(h, d and lambda)
## this can be fitted on a laptop: default number of thread = 1

## example 1: estimated graph at 265th time point (2008-01-22) using pre-tuned hyperparameters (using method "loggle")
## time cost: 1062.77s
ts <- proc.time()
ind <- 265
pos.ind <- which(pos == ind)
result.1 <- loggle.cv.vote(X, ind, h = tuned.parameter$loggle$h, d = tuned.parameter$loggle$d[pos.ind], 
                           lambda = tuned.parameter$loggle$lambda[pos.ind], fit.type = "pseudo", epi.abs = 1e-4, 
                           epi.rel = 1e-2, num.thread = 1)
te <- proc.time()
sprintf("Time used for model fitting at %dth time point: %.2fs", ind, (te-ts)[3])
sprintf("Number of edges at %dth time point: %d", ind, result.1$edge.num)
# Number of edges at 265th time point: 511

## example 2: estimated graph at 743rd time point (2009-12-11) using pre-tuned hyperparameters (using method "loggle")
## time cost: 998.07s
ts <- proc.time()
ind <- 743
pos.ind <- which(pos == ind)
result.2 <- loggle.cv.vote(X, ind, h = tuned.parameter$loggle$h, d = tuned.parameter$loggle$d[pos.ind], 
                           lambda = tuned.parameter$loggle$lambda[pos.ind], fit.type = "pseudo", epi.abs = 1e-4, 
                           epi.rel = 1e-2, num.thread = 1)
te <- proc.time()
sprintf("Time used for model fitting at %dth time point: %.2fs", ind, (te-ts)[3])
sprintf("Number of edges at %dth time point: %d", ind, result.2$edge.num)
# Number of edges at 743th time point: 841