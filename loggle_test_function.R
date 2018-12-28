
# Evaluation function for loggle #######################################################################################
########################################################################################################################

# Input ###
# pos: indices of time points where graphs are estimated
# Omega: a list of precision matrices or adjacency matrices across time points specified by pos
# Omega.true: a list of true precision matrices across all time points
# X: a p by N data matrix containing observations on a time grid ranging from 0 to 1

# Output ###
# FDR: mean and standard deviation of false discovery rates across time points specified by pos
# power: mean and standard deviation of powers across time points specified by pos
# F1: mean of F1 scores across time points specified by pos
# KLdis: mean and standard deviation of KL distance across time points specified by pos

loggle.eval <- function(pos, Omega, Omega.true, X) {
  
  K <- length(pos)
  p <- dim(X)[1]
  tri.ind <- lower.tri(Omega.true[[1]])
  
  edge.true <- rep(NA, K)
  edge <- rep(NA, K)
  overlap.edge <- rep(NA, K)
  
  for(k in 1:K) {
    adj.mat.true <- as.matrix(Omega.true[[pos[k]]] != 0)
    adj.mat <- as.matrix(Omega[[k]] != 0)
    edge.true[k] <- sum(adj.mat.true[tri.ind])
    edge[k] <- sum(adj.mat[tri.ind])
    overlap.edge[k] <- sum((adj.mat.true&adj.mat)[tri.ind])
  }
  
  FDR.list <- 1-overlap.edge/edge
  FDR.list[is.nan(FDR.list)] <- 0
  power.list <- overlap.edge/edge.true
  
  FDR <- c(mean(FDR.list), sd(FDR.list))
  power <- c(mean(power.list), sd(power.list))
  F1 <- mean(2*(1-FDR.list)*power.list/(1-FDR.list+power.list))
  
  KLdis.list <- rep(NA, K)
  Omega.refit <- loggle.refit(X, pos, Omega, print.detail = FALSE)
  for(k in 1:K) {
    Omega.prod <- t(solve(as.matrix(Omega.true[[pos[k]]]), as.matrix(Omega.refit[[k]])))
    KLdis.list[k] <- sum(diag(Omega.prod)) - log(det(Omega.prod)) - p
  }
  KLdis <- c(mean(KLdis.list), sd(KLdis.list))
  
  result <- list(FDR = FDR, power = power, F1 = F1, KLdis = KLdis)
  return(result)
}


# Selection function for cross validation result for method "loggle", "kernel" and "invariant" #########################
########################################################################################################################

# Input ###
# cv.result: results from loggle.cv
# select.type: "all_flexible": optimal d and lambda can vary across time points specified by pos,
#              "d_fixed": optimal d is fixed and optimal lambda can vary across time points specified by pos,
#              "all_fixed": optimal d and lambda are fixed across time points specified by pos
# cv.vote.thres: an edge is kept after cv.vote if and only if it exists in no less than cv.vote.thres*cv.fold cv folds
# method: "loggle": loggle method, 
#         "kernel": Zhou's method,
#         "invariant": Wang's method

# Output ###
# h.opt: optimal value of h
# d.opt: a vector of optimal values of d for each estimated graph
# lambda.opt: a vector of optimal values of lambda for each estimated graph
# cv.score.opt: optimal cv score (averaged over time points and cv folds)
# edge.num.opt: a vector of numbers of edges for each estimated graph
# edge.opt: a list of edges for each estimated graph
# adj.mat.opt: a list of adjacency matrices for each estimated graph

loggle.cv.select.full <- function(cv.result, select.type = "all_flexible", cv.vote.thres = 0.8, method = "loggle") {
  
  d.list <- as.numeric(colnames(cv.result$cv.result.h[[1]]$cv.score))
  H <- length(cv.result$cv.result.h)
  cv.fold <- length(cv.result$cv.result.h[[1]]$cv.result.fold)
  
  if(method %in% c("kernel", "invariant")) {
    d.pos <- ifelse(method == "kernel", which(d.list == 0), which(d.list == 1))
    for(h in 1:H) {
      cv.result$cv.result.h[[h]]$cv.score <- cv.result$cv.result.h[[h]]$cv.score[ , d.pos, , , drop=F]
      for(i in 1:cv.fold) {
        cv.result$cv.result.h[[h]]$cv.result.fold[[i]]$Omega <- 
          cv.result$cv.result.h[[h]]$cv.result.fold[[i]]$Omega[ , d.pos, , drop=F]
      }
    }
  } else if(method != "loggle") {
    stop("method must be 'loggle', 'kernel' or 'invariant'!")
  }
  
  cv.select.result <- loggle.cv.select(cv.result, select.type, cv.vote.thres)
  return(cv.select.result)
}


# Summary matrix containing number of within-sector and cross-sector edges #############################################
########################################################################################################################

# Input ###
# select.result: result from loggle.cv.select.full
# sp.num: number of stocks within each sector
# h: bandwidth in kernel smoothed edge number

# Output ###
# edge.sector: row 1 to row length(sp.num): within-sector edge number at each time point for each sector
#              row length(sp.num)+1: cross-sector edge number at each time point
#              row length(sp.num)+2: total edge number at each time point
# edge.sector.sm: smoothed edge.sector across time points with kernel bandwidth h

loggle.sector <- function(select.result, sp.num, h) {
  
  K <- length(select.result$edge.num.opt)
  sector.num <- length(sp.num)
  sp.ind <- c(0, cumsum(sp.num))
  edge.sector <- matrix(NA, sector.num+2, K)
  edge.sector.sm <- matrix(NA, sector.num+2, K)
  
  for(i in 1:sector.num) {
    edge.sector[i, ] <- sapply(1:K, function(k) 
      (sum(select.result$adj.mat.opt[[k]][(sp.ind[i]+1):sp.ind[i+1],(sp.ind[i]+1):sp.ind[i+1]])-sp.num[i])/2)
  }
  edge.sector[sector.num+2, ] <- select.result$edge.num.opt
  edge.sector[sector.num+1, ] <- edge.sector[sector.num+2, ] - colSums(edge.sector[1:sector.num, ])
  
  for(k in 1:K) {
    Kh <- pmax(3/4*(1-((1:K-k)/((K-1)*h))^2), 0)
    omega <- Kh/sum(Kh)
    edge.sector.sm[, k] <- edge.sector %*% omega
  }

  result <- list(edge.sector = edge.sector, edge.sector.sm = edge.sector.sm)
  return(result)
}