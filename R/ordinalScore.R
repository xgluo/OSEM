library(dplyr)
library(leaps)
library(pcalg)
library(glmnet)

##' learnCuts(param):
##' a function that estimates the thresholds/cuts marginally for each variable
##' @param param: an ordinal scoreparameters object
##' @return a list of cuts of length n
learnCuts <- function(param) {
  cuts <- list()
  for (i in c(1:param$n)) {
    tb <- table(param$data[,i])
    tblevels <- names(tb)
    temptb <- c()
    if (length(tblevels) <- param$ordinalLevels[i]) {
      for (a in as.character(c(1:param$ordinalLevels[i])-1)) {
        if (!(a %in% tblevels)) {
          temptb[a] <- 1
        } else {
          temptb[a] <- tb[a] + 1
        }
      }
      tb <- cumsum(temptb) / sum(temptb)
    } else {
      tb <- cumsum(tb) / param$N
    }
    cuts <- append(cuts, list(c(-Inf, sapply(tb, function (x) qnorm(x)))))
  }
  return(cuts)
}


##' getPairwiseTb(param,i,j):
##' a function that gets the contingency table of two ordinal variables with pre-specified number of levels
##' @param param: an ordinal scoreparameters object
##' @param i: variable X_i
##' @param j: variable X_j
##' @return estimated table between X_i and X_j
getPairwiseTb <- function(param, i, j) {

  tb <- matrix(1, nrow = param$ordinalLevels[i], ncol = param$ordinalLevels[j])

  for (k in c(1:param$N)) {
    x <- param$data[k,i]
    y <- param$data[k,j]
    tb[x+1, y+1] <- tb[x+1, y+1] + 1
  }

  dimnames(tb) <- list(as.character(c(1:param$ordinalLevels[i]) - 1),
                       as.character(c(1:param$ordinalLevels[j]) - 1))
  return(tb)

}


##' learnRhoij(param,i,j):
##' a function that estimates the latent Gaussian correlation between two binary variables
##' @param param: an ordinal scoreparameters object
##' @param i: variable X_i
##' @param j: variable X_j
##' @return estimated correlation between X_i and X_j
learnRhoij <- function(param,i,j) {

  tb <- getPairwiseTb(param, i, j)
  cutsfori <- simplify2array(param$cuts[i])
  cutsforj <- simplify2array(param$cuts[j])
  Li <- param$ordinalLevels[i]
  Lj <- param$ordinalLevels[j]

  # an objective function that computes the log-likelihood
  objective <- function(rho) {
    LL <- 0
    for (li in c(1:Li)) {
      for (lj in c(1:Lj)) {
        if (tb[li,lj] != 0) {
          prob <- mvtnorm::pmvnorm(lower = c(cutsfori[li], cutsforj[lj]),
                                   upper = c(cutsfori[li+1], cutsforj[lj+1]),
                                   corr = matrix(c(1,rho,rho,1),nrow = 2))[1]
          LL <- LL + log(prob) * tb[li,lj]
        }
      }
    }
    return(LL)
  }

  # search for rho that maximizes the pairwise log-likelihood
  rho_hat <- optimize(objective, interval = c(-1 + 1e-3,1 - 1e-3),maximum = TRUE)
  return(rho_hat$maximum)
}

##' refine_parameters(param):
##' a helper function that refines the covariance matrix based on observed data
##' @param param: an ordinal scoreparameters object
##' @return an updated ordinal scoreparameters object
refineParameters <- function(param) {

  for (i in c(1:5)) {
    old_Sigma_hat <- param$Sigma_hat
    param$hidden_data <- sampleTMVN(param)
    param$Sigma_hat <- cor(param$hidden_data)
    param_diff <- sum((param$Sigma_hat - old_Sigma_hat)^2)
    if (param_diff < 0.01) {
      return(param)
    } else {
      cat("diff in param: ", param_diff, "\n")
    }
  }
  return(param)

}


##' ordinalScoreParam(initparam,ordType): a function to be called by 'scoreparameters' in BiDAG
##' It prepares the parameters needed to calculate the score for ordinal data
##' @param initparam: parameters returned from 'scoreparameters'
##' @param penType: type of penalization to be applied ("AIC" or "BIC" or "other")
##' @param L: number of truncated multivariate copies to be generated for each observation
##' @param lambda: coefficient for the BIC penalty term (default: 0 no penalty)
##' @return the input parameters
##' + estimated thresholds
##' + estimated latent Gaussian correlation matrix
ordinalScoreParam <- function(initparam,
                              usrpar = list(penType = c("AIC","BIC","other"),
                                            L = 5,
                                            lambda = 2,
                                            preLevels = NULL)) {

  n <- initparam$n
  initparam$penType <- usrpar$penType
  initparam$lambda <- usrpar$lambda
  initparam$L <- usrpar$L
  initparam$N <- nrow(initparam$data)

  # Convert to matrix/array if necessary
  if (!is.matrix(initparam$data)) {
    if (is.data.frame(initparam$data)) {
      initparam$data <- as.matrix(initparam$data)
    } else {
      stop("Ordinal data table is not a matrix")
    }
  }

  # Determine ordinal levels
  if (is.null(usrpar$preLevels)) {
    # Reassign levels if necessary
    for (i in c(1:n)) {
      x <- sort(unique(initparam$data[,i]))
      if (!setequal(x, c(0:max(x)))) {
        cat("Some levels in variable",i,"are not present in the data, merging levels...\n")
        for (j in c(1:length(x))) {
          ind <- which(initparam$data[,i] == x[j])
          initparam$data[ind,i] <- j - 1
        }
      }
    }

    # Count the number of levels for each variable
    initparam$ordinalLevels <- apply(initparam$data, 2, function (x) length(unique(x)))
  } else {
    initparam$ordinalLevels <- usrpar$preLevels
  }

  if (min(initparam$ordinalLevels) < 2) {
    stop("All variables need to have at least two levels...")
  }

  initparam$hidden_data <- initparam$data

  # Learn the cuts/thresholds
  initparam$cuts <- learnCuts(initparam)

  # Estimate Gaussian correlation matrix
  Sigma_hat <- diag(initparam$n)

  cat("Estimating correlation matrix...\n")
  for (i in c(1:(n-1))) {
    for (j in c((i+1):n)) {
      Sigma_hat[i,j] <- learnRhoij(initparam,i,j)
    }
  }
  Sigma_hat <- Sigma_hat + t(Sigma_hat * upper.tri(Sigma_hat))
  dimnames(Sigma_hat) <- list(c(1:n),c(1:n))

  cat("Smallest eigenvalue of the estimated correlation matrix is", min(eigen(Sigma_hat)$values), "...\n")
  # Smooth the correlation matrix to make it positive definite
  # Usually necessary for large DAGs because of pairwise estimation
  # Only done when the smallest eigenvalue is negative
  initparam$Sigma_hat <- psych::cor.smooth(Sigma_hat,eig.tol = 1e-3)

  print("Refining parameters...")
  initparam <- refineParameters(initparam)

  return(initparam)
}

##' ordinalCoreScore(j,parentnodes,n,param): a function to be called by 'DAGcorescore' in BiDAG
##' It computes the ordinal score of a node given its parents
##' @param j: index of node j
##' @param parentnodes: indices of parent nodes of j
##' @param n: dimension of the variable
##' @param param: parameters returned from 'scoreparameters'
##' @return the ordinal score of a node given its parents
ordinalCoreScore <- function(j,parentnodes,n,param) {

  if (j %in% parentnodes) {
    stop("The parent set contains the current node.")
  }

  lp <- length(parentnodes) # number of parents
  corescore <- 0
  N <- param$N

  switch(as.character(lp),
         "0" = { # no parents
           corescore <- -N/2 * (1 + log(2 * pi))
         },
         "1" = { # one parent
           corescore <- -N/2 * (log(param$Sigma_hat[j,j] - param$Sigma_hat[j,parentnodes]^2 / param$Sigma_hat[parentnodes,parentnodes]) + 1 + log(2 * pi))
         },
         { # more parents
           b <- param$Sigma_hat[j,parentnodes]
           S_p <- param$Sigma_hat[parentnodes,parentnodes]
           choltemp <- chol(S_p)
           corescore <- -N/2 * (log(param$Sigma_hat[j,j] - sum(backsolve(choltemp,b,transpose=TRUE)^2)) + 1 + log(2 * pi))
         }
  )

  if (param$penType == "AIC") {
    param$lambda <- 2 / log(param$N)
  } else if (param$penType == "BIC") {
    param$lambda <- 1
  }

  return(corescore - param$lambda * log(param$N) / 2 * lp)
}

##' sampleTMVN(param):
##' a function that samples truncated multivariate Gaussian data
##' @param param: scoreparameters object
##' @return N times L number of truncated multivariate Gaussian data points
sampleTMVN <- function(param) {
  Y <- NULL
  # Obtain different combinations and counts for ordinal data
  comb_X <- suppressMessages(data.frame(param$data) %>% dplyr::group_by_all() %>% dplyr::summarise(COUNT = n()))
  for (i in (1:nrow(comb_X))) {
    x <- as.numeric(comb_X[i,(1:param$n)])
    lb <- sapply(c(1:param$n), function (j) simplify2array(param$cuts[j])[x[j] + 1])
    ub <- sapply(c(1:param$n), function (j) simplify2array(param$cuts[j])[x[j] + 2])
    init <- sapply(c(1:param$n), function (j) tmvmixnorm::rtuvn(lower = lb[j], upper = ub[j]))
    Y <- rbind(Y, tmvmixnorm::rtmvn(n = param$L * as.numeric(comb_X[i,(param$n+1)]),
                                    Mean = rep(0, param$n),
                                    Sigma = param$Sigma_hat,
                                    lower = lb,
                                    upper = ub,
                                    int = init))
  }
  return(Y)
}

##' getExpectedStats(param):
##' a function that takes the score parameters, samples latent Gaussian data
##' and computes expected statistics (mean and correlation matrix)
##' @param param: a scoreparameter object
##' @return an updated score parameters with the sampled latent data and expected statistics
getExpectedStats <- function(param) {

  print("Updating expected statistics...")
  param$hidden_data <- sampleTMVN(param)
  param$Sigma_hat <- cor(param$hidden_data)
  return(param)

}

##' ordinalUpdateParam(param,AM):
##' a function that estimates the Gaussian mean and correlation matrix given a sample of DAGs
##' @param param: a scoreparameter object for ordinal data
##' @param AM: an adjacency matrix of a DAG
##' @param regsubsets: a boolean value indicating whether subset selection is used in parameter update (default: TRUE)
##' @return an ordinal scoreparameter object with updated Gaussian correlation matrix
ordinalUpdateParam <- function(param,AM,regsubsets = TRUE) {

  start_time = Sys.time()

  cat("Estimating parameters...\n")

  if (!is.matrix(AM)) {
    if (class(AM)  == "graphNEL") {
      AM <- as(AM,"matrix")
    } else {
      stop("The parameter AM needs to be an adjacency matrix...")
    }
  }

  if (nrow(AM) != param$n) {
    stop("The adjacency matrix provided has wrong dimension...")
  }

  Y <- param$hidden_data
  n <- param$n
  N <- param$N
  NL <- N * param$L
  sorted_ind <- ggm::topOrder(AM)
  var_temp <- rep(1,n)
  B <- matrix(0,nrow = n, ncol = n)

  # Estimate the conditional parameters following the topological order
  for (i in (1:n)) {

    j <- sorted_ind[i]
    parentnodes <- which(AM[,j] == 1)
    lp <- length(parentnodes)
    Y_j <- Y[,j]

    switch(as.character(lp),

           "0" = { # no parents

             var_temp[j] <- mean(Y_j^2)

           },

           "1" = { # one parent

             X <- matrix(Y[,parentnodes])
             var_temp[j] <- mean(Y_j^2)
             add1 <- lm(Y_j ~ X - 1)
             add1MSE <- mean(add1$residuals^2)

             if (regsubsets) {

               if (NL*log(add1MSE) + param$lambda * log(NL) < NL*log(var_temp[j])) {
                 beta <- add1$coef
                 B[j,parentnodes] <- beta
                 var_temp[j] <- add1MSE
               }

             } else {
               beta <- add1$coef
               B[j,parentnodes] <- beta
               var_temp[j] <- add1MSE
             }

           },

           { # more parents

             X <- Y[,parentnodes]

             if (regsubsets) {

               regfit.full <- regsubsets(x = X, y = Y_j,nvmax = lp,nbest = 1,
                                         method = "exhaustive",intercept = FALSE)
               reg.summary <- summary(regfit.full)
               reg.bic <- reg.summary$bic

               if (param$penType == "other") {
                 reg.bic <- reg.bic - (1 - param$lambda) * log(NL) * apply(reg.summary$which,1,sum)
               }

               best.n <- which.min(reg.bic)

               beta <- coef(regfit.full, id = best.n)
               B[j,parentnodes[reg.summary$which[best.n,]]] <- beta
               var_temp[j] <- reg.summary$rss[best.n] / NL

             } else {

               lm.full <- lm(Y_j ~ X - 1)
               B[j,parentnodes] <- coefficients(lm.full)
               var_temp[j] <- sum(resid(lm.full)^2) / NL

             }

           })

  }

  I_B <- diag(n) - B
  I_B.inv <- solve(I_B)
  V <- diag(var_temp)
  param$Sigma_hat <- cov2cor(I_B.inv %*% V %*% t(I_B.inv))

  end_time = Sys.time()
  print(end_time - start_time)

  return(param)
}

######### Taken from the source code of the pcalg package
## this function takes as parameter the adjacency matrix of a pdag (or cpdag)
## and returns the pattern of this pdag in the Meek sense, that is,
## it returns the adjacency matrix of the graph with the same skeleton where the only oriented
## edges are the v-structures (can be easily modified to work for MAGs/PAGs)
getPattern <- function(amat){

  ## makes the whole graph undirected
  tmp <- amat + t(amat)
  tmp[tmp == 2] <- 1

  ## find all v-structures i -> k <- j s.t. i not adj to k
  ## and make only those edges directed
  for (i in 1: (length(tmp[1,])-1)) {
    for (j in (i+1): length(tmp[1,])){
      if ((amat[j,i] ==0) & (amat[i,j] ==0) & (i!=j)){ ## if i no adjacent with j in G

        possible.k <- which(amat[i,]!= 0 & amat[,i]==0) ## finds all k such that i -> k is in G

        if (length(possible.k)!=0){    ## if there are any such k's then check whether j -> k for any of them
          for (k in 1: length(possible.k)){
            if ((amat[j,possible.k[k]] ==1) & (amat[possible.k[k],j]==0)) { ## if j -> k add the v-struc orientation to tmp
              tmp[possible.k[k],i] <- 0
              tmp[possible.k[k],j] <- 0
            }
          }
        }
      }
    }
  }
  tmp
}

convert_to_skeleton <- function(DAG) {
  DAG <- DAG + t(DAG)
  return(DAG != 0)
}


##' comparePatterns(estDAG,trueDAG):
##' a function that compares the patterns of two DAGs in the sense of Meek (1995)
##' @param estDAG: estimated DAG (need to be a matrix or a graphNEL object)
##' @param trueDAG: true DAG (need to be a matrix or a graphNEL object)
##' @param hardP2P: (default: FALSE) An edge in the estimated pattern is counted as 1 true positive,
##' if it has exactly the same direction (directed/undirected) as the corresponding edge in the true pattern.
##' Otherwise, it is counted as 1 false positive. When FALSE, an edge in the estimated pattern
##' is counted as 0.5 true positive and 0.5 false positive,
##' if exactly one of this edge and the corresponding edge in the true pattern is undirected.
##' @return an array of metrics
comparePatterns <- function(estDAG, trueDAG, hardP2P = FALSE, skeleton = FALSE) {

  if (skeleton) { # use skeletons
    trueSkel <- convert_to_skeleton(as(trueDAG, "matrix"))
    estSkel <- convert_to_skeleton(as(estDAG, "matrix"))

    temp1 <- estSkel[upper.tri(estSkel)]
    temp2 <- trueSkel[upper.tri(trueSkel)]

    pred_P <- sum(temp1 != 0)
    true_P <- sum(temp2 != 0)
    true_N <- sum(temp2 == 0)

    TP <- sum((temp1 != 0) * (temp2 != 0))
    FP <- pred_P - TP
    FN <- sum((temp1 == 0) * (temp2 != 0))
    TN <- sum((temp1 == 0) * (temp2 == 0))

  } else { # use patterns
    # Convert estimated DAG to CPDAG
    if (is.matrix(estDAG)) {
      estPDAG <- dag2cpdag(estDAG)
    } else if (class(estDAG) == "pcAlgo") {
      estPDAG <- t(as(estDAG,"matrix"))
    } else if (class(estDAG) == "graphNEL") {
      estPDAG <- dag2cpdag(as(estDAG,"matrix"))
    }

    # Convert true DAG to CPDAG
    if (is.matrix(trueDAG)) {
      truePDAG <- dag2cpdag(trueDAG)
    } else if (class(trueDAG) == "graphNEL") {
      truePDAG <- dag2cpdag(as(trueDAG,"matrix"))
    } else {
      stop("Please check if the true DAG is indeed a DAG...")
    }

    # Convert CPDAGs to patterns
    truePattern <- getPattern(truePDAG)
    estPattern <- getPattern(estPDAG)

    # 0: no edge; 1: -->; 2: <--; 3: <-->
    temp1 <- estPattern[upper.tri(estPattern)] + 2 * t(estPattern)[upper.tri(t(estPattern))]
    temp2 <- truePattern[upper.tri(truePattern)] + 2 * t(truePattern)[upper.tri(t(truePattern))]

    # Number of edges in the estimated pattern
    pred_P <- sum(temp1 != 0)

    # Number of edges in the true pattern
    true_P <- sum(temp2 != 0)

    # Number of non-edges in the true pattern
    true_N <- sum(temp2 == 0)

    # TP, FP, TN, FN, SHD
    if (hardP2P) {
      TP <- sum((temp1 != 0) * (temp1 == temp2))
    } else {
      TP <- sum((temp1 != 0) * (temp1 == temp2)) + 0.5 * sum((temp1 * temp2) == 3) + 0.5 * sum((temp1 * temp2) == 6)
    }
    FP <- pred_P - TP
    FN <- sum((temp1 == 0) * (temp2 != 0))
    TN <- sum((temp1 == 0) * (temp2 == 0))
  }

  # Structural hamming distance
  SHD <- FP + FN

  # Precision
  if ((TP + FP) == 0) {
    Precision <- 0
  } else {
    Precision <- TP / (TP + FP)
  }

  # TPR, FPR_P, FPR_N
  if (true_P == 0) { # true graph is empty
    if (FP >= 0) {
      TPR <- 0
      FPR_P <- 1
    } else {
      TPR <- 1
      FPR_P <- 0
    }
  } else { # true graph is non-empty
    TPR <- TP / true_P
    FPR_P <- FP / true_P
  }

  if (true_N == 0) { # true graph is full
    FPR_N <- 0
  } else { # true graph is not full
    FPR_N <- FP / true_N
  }

  compPattern <- c(SHD,TP,FP,TN,FN,Precision,TPR,FPR_N,FPR_P)
  names(compPattern) <- c("SHD","TP","FP","TN","FN","Precision","TPR","FPR_N","FPR_P")
  return(round(compPattern,2))
}

##' getSHD(estDAG, trueDAG):
##' a function that computes the structural difference between the PDAGs of two DAGs
##' @param estDAG: adjacency matrix of the estimated DAG
##' @param trueDAG: adjacency matrix of the true DAG
##' @return structural difference between two PDAGs
getSHD <- function(estDAG, trueDAG) {

  if (!is.matrix(estDAG)) {
    estDAG <- as(estDAG,"matrix") * 1
  }

  if (!is.matrix(trueDAG)) {
    trueDAG <- as(trueDAG,"matrix") * 1
  }

  DAGdiff <- dag2cpdag(estDAG) != dag2cpdag(trueDAG)
  return(sum(as.logical(DAGdiff + t(DAGdiff)))/2)
}

##' observedLL(param):
##' a function that computes observed log-likelihood in the ordinal setting
##' @param param: a scoreparameter object for ordinal data
##' @return observed log-likelihood
observedLL <- function(param) {
  LL <- 0
  comb_X <- suppressMessages(data.frame(param$data) %>% dplyr::group_by_all() %>% dplyr::summarise(COUNT = n()))
  for (i in (1:nrow(comb_X))) {
    x <- as.numeric(comb_X[i,(1:param$n)])
    lb <- sapply(c(1:param$n), function (j) simplify2array(param$cuts[j])[x[j] + 1])
    ub <- sapply(c(1:param$n), function (j) simplify2array(param$cuts[j])[x[j] + 2])
    prob <- mvtnorm::pmvnorm(lower = lb,
                             upper = ub,
                             corr = param$Sigma_hat)[1]
    LL <- LL + as.numeric(comb_X[i,(param$n+1)]) * log(prob)
  }
  return(LL)
}

##' ordinalStructEM(n, data, usrpar):
##' a function that implements the structural EM algorithm for learning
##' Bayesian networks from ordinal data
##' @param n: dimension of the variable
##' @param data: ordinal dataset
##' @param usrpar: list of parameters for ordinal scores
##' @param computeObservedLL: compute the observed likelihood of the parameters (default: FALSE)
##' @param iterMCMC_alpha: significance level for the constraint-based initialization (default: 0.05)
##' @param regsubsets: a boolean value indicating whether subset selection is used in parameter update (default: TRUE)
##' @return an object of class MCMCtrace
ordinalStructEM <- function(n, data,
                            usrpar = list(penType = c("AIC","BIC","other"),
                                          L = 5,
                                          lambda = 2,
                                          preLevels = NULL),
                            computeObservedLL = FALSE,
                            iterMCMC_alpha = 0.05,
                            regsubsets = TRUE) {

  start_time = Sys.time()
  print("Initializing parameters and DAG...")
  param <- scoreparameters("usr",data,usrpar = usrpar)

  print("SEM iteration starts...")
  SHD <- 1
  iter <- 0
  currentBestDAG <- NULL

  # Maximum iterations to control runtime (can change)
  if ((param$lambda <= 1) || (n >= 30)) {
    nr_EM_iter <- 5
    nr_plus1it <- 5
  } else {
    nr_EM_iter <- 10
    nr_plus1it <- 10
  }

  while ((SHD != 0) && (iter < nr_EM_iter)) {

    # Compute expected statistics
    param <- getExpectedStats(param)

    # Structure update
    currentDAGobj <- iterativeMCMC(param, plus1it = nr_plus1it, hardlimit = n, alpha = iterMCMC_alpha)
    candidateBestDAG <- as.matrix(currentDAGobj$DAG)

    # Parameter update
    param <- ordinalUpdateParam(param,candidateBestDAG, regsubsets = regsubsets)

    # Check convergence
    if (!(is.null(currentBestDAG))) {
      SHD <- getSHD(candidateBestDAG, currentBestDAG)
      cat("PDAG SHD: ", SHD, "\n")
    }
    currentBestDAG <- candidateBestDAG

    iter <- iter + 1
  }
  print("SEM iteration ends...")

  if (computeObservedLL) {
    print("Calculating observed log-likelihood...")
    currentDAGobj$observed_score <- observedLL(param) - param$lambda * log(param$N) / 2 * sum(currentBestDAG)
  }

  currentDAGobj$param <- param

  end_time = Sys.time()
  print(end_time - start_time)
  currentDAGobj$runtime <- as.double(end_time - start_time,units = "secs")

  return(currentDAGobj)

}

##' BGe_MAP_Sigma(R, DAG):
##' a function that computes the MAP Covariance matrix based on a DAG estimated
##' using the BGe score
##' Source: https://arxiv.org/pdf/2010.00684.pdf (Equation 2)
##' @param R: BGe MAP covariance matrix
##' @param DAG: a DAG estimated using the BGe score
##' @return A DAG-aware BGe MAP covariance matrix
BGe_MAP_Sigma <- function(R, DAG) {

  n <- nrow(DAG)
  B <- matrix(0, nrow = n, ncol = n)

  for (j in c(1:ncol(DAG))) {
    pa <- which(DAG[,j] == 1)
    lp <- length(pa)
    switch (as.character(lp),
            "0" = {
              next
            },
            "1" = {
              B[j,pa] <- 1 / R[pa,pa] * R[pa,j]
            },
            {
              B[j,pa] <- chol2inv(chol(R[pa,pa])) %*% R[pa,j]
            }
    )
  }

  I_B <- diag(n) - B
  I_B.inv <- solve(I_B)
  V <- diag(diag(R))
  Sigma <- cov2cor(I_B.inv %*% V %*% t(I_B.inv))

  return(Sigma)
}

##' rmvDAG2(N, randDAGobj):
##' a function that does the same thing as the pcalg::rmvDAG function
##' but the input DAG is not necessarily topologically ordered
##' @param N: number of samples to be drawn
##' @param randDAGobj: a graph object generated from the pcalg::randDAG function
##' @return a Gaussian dataset
rmvDAG2 <- function(N, randDAGobj) {
  AM <- as(randDAGobj, "matrix")
  sorted_ind <- ggm::topOrder((AM != 0))
  n <- nrow(AM)
  data <- matrix(nrow = N,ncol = n)
  for (j in sorted_ind) {
    parentnodes <- which(AM[,j] != 0)
    lp <- length(parentnodes)
    switch (as.character(lp),
            "0" = {data[,j] <- rnorm(N)},
            "1" = {data[,j] <- rnorm(N, mean = data[,parentnodes] * AM[parentnodes,j], sd = 1)},
            {data[,j] <- rnorm(N, mean = data[,parentnodes] %*% AM[parentnodes,j], sd = 1)}
    )
  }
  return(data)
}

##' cutfun(L, c):
##' a function that simulates the cell probabilities from a symmetric Dirichlet distribution
##' @param L: number of ordinal levels
##' @param c: Dirichlet concentration parameter
##' @return a list of probabilities of length L, summing up to 1
cutfun <- function(L,c) {
  p <- gtools::rdirichlet(1,rep(c,L))
  return(qnorm(cumsum(p)[1:(L-1)]))
}

##' convertToOrdinal(scaled_data, exp_levels, concent_param):
##' a function that converts standardized Gaussian data into ordinal data
##' @param scaled_data: Gaussian dataset with each dimension standardized
##' @param exp_levels: expected number of ordinal levels
##' @param concent_param: Dirichlet concentration parameter
##' @return an ordinal dataset
convertToOrdinal <- function(scaled_data, exp_levels = 4,concent_param = 2) {
  n <- ncol(scaled_data)
  if (exp_levels == 2) {
    ordinal_levels <- replicate(n,2)
  } else {
    ordinal_levels <- replicate(n,sample(c(2:(2 * exp_levels - 2)),1))
  }
  ordinal_data <- scaled_data
  for (i in c(1:n)) {

    check_levels <- ordinal_levels[i] - 1
    while (check_levels != ordinal_levels[i]) {
      cuts <- c(-Inf,
                cutfun(ordinal_levels[i],concent_param),
                Inf)
      temp <- cut(scaled_data[,i], simplify2array(cuts), labels = FALSE) - 1
      check_levels <- length(unique(temp))
    }
    ordinal_data[,i] <- temp

  }
  colnames(ordinal_data) <- c(1:n)
  return(ordinal_data)
}

##' mywFUN(m):
##' a function that samples the edge weights uniformly from the interval (-1,-0.4) U (0.4,1)
##' @param m: number of edges in the DAG
##' @return m edge weights
mywFUN <- function(m) {
  return(replicate(m,mywFUNhelper()))
}
mywFUNhelper <- function() {
  y <- runif(1, 0, 1.2)
  if( y < 0.6 ){
    x <- -1 + y
  }else{
    x <- 0.4 + y - 0.6
  }
  return(x)
}

generateOrdinal <- function(N, n, trueDAG, exp_levels = 4, concent_param = 2) {
  
  hidden_data <- rmvDAG2(N, trueDAG)
  scaled_data <- t(t(hidden_data) - apply(hidden_data,2,mean))
  truecov <- trueCov(trueDAG)
  D <- diag(sqrt(diag(truecov)))
  D.inv <- chol2inv(chol(D))
  trueSigma <- D.inv %*% truecov %*% D.inv
  scaled_data <- t(D.inv %*% t(scaled_data))
  
  # Convert the Gaussian dataset into an ordinal dataset
  ordinal_data <- convertToOrdinal(scaled_data, exp_levels = exp_levels, concent_param = concent_param)
  
  return(ordinal_data)
  
}
