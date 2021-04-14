### 1. Load necessary packages and source code

``` r
library(BiDAG)
library(Rgraphviz)
library(pcalg)
library(graph)
library(MXM)
library(sbgcop)
library(infotheo)
source('~/Documents/Projects/OSEM/CausalMissingValues/R/gaussCItestLocal.R')
source('~/Documents/Projects/OSEM/CausalMissingValues/R/inferCopulaModel.R')
source('~/Documents/Projects/OSEM/CausalMissingValues/R/addMAR.R')
```

``` r
setwd("../R")
# Major file containing the OSEM algorithm
source("ordinalScore.R") 
# Modify some of the existing functions in the BiDAG package to accommodate our user-defined functions
insertSource("spacefns.R",package = "BiDAG")
insertSource("usrscorefns.R",package = "BiDAG")
```

### 2. Other functions

``` r
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
```

``` r
##' cutfun(L, c): 
##' a function that simulates the cell probabilities from a symmetric Dirichlet distribution
##' @param L: number of ordinal levels
##' @param c: Dirichlet concentration parameter
##' @return a list of probabilities of length L, summing up to 1
cutfun <- function(L,c) {
  p <- gtools::rdirichlet(1,rep(c,L))
  return(qnorm(cumsum(p)[1:(L-1)]))
}
```

``` r
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
```

``` r
##' catCItest(x, y, S, suffStat): 
##' a helper function for categorical tests to be consumed by the pcalg::pc function
##' @param x,y,S: position of variable X, Y and set of variables S, respectively, in the adjacency matrix
##' @param suffStat: sufficient statistics for the test
##' @return p-value of the categorical tests
catCItest <- function(x, y, S, suffStat) {
  d <- as.data.frame(suffStat$data)
  d[] <- lapply(d[], as.ordered)
  return(gRim::ciTest_ordinal(d,
                              set = as.numeric(c(x,y,S)),
                              statistic = suffStat$stat_type)$P)
}
```

``` r
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
```

### 3. Generate DAGs and ordinal data

``` r
set.seed(888)
# Generate a regular DAG with 20 nodes with 4 number of neighbors
n <- 20
trueDAG <- randDAG(n = n, d = 4, method = "power", wFUN = list(mywFUN))

# Generate a Gaussian dataset based on the DAG and standardize each dimension
# Sample size = 500
N <- 500
hidden_data <- rmvDAG2(N, trueDAG)
truecov <- trueCov(trueDAG)
D <- diag(sqrt(diag(truecov)))
D.inv <- chol2inv(chol(D))
trueSigma <- D.inv %*% truecov %*% D.inv
scaled_data <- t(D.inv %*% t(hidden_data))

# Convert the Gaussian dataset into an ordinal dataset
ordinal_data <- convertToOrdinal(scaled_data, exp_levels = 4,concent_param = 2)
```

### 4. Learn the structures using different methods

The following simulations can be repeated many times for different
configurations of parameters in order to obtain the ROC curves.

``` r
# NPC algorithm with the nominal deviance test (significance level: 0.05)
NPCfit <- pc(suffStat = list(data = ordinal_data, stat_type = "dev"),
             alpha = 0.05,
             indepTest = catCItest,
             labels = colnames(ordinal_data))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# Compare the patterns between them
comparePatterns(NPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     22.00      7.00      7.00    161.00     15.00      0.50      0.28      0.04 
    ##     FPR_P 
    ##      0.28

``` r
comparePatterns(NPCfit,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     20.50      8.50      5.50    161.00     15.00      0.61      0.34      0.03 
    ##     FPR_P 
    ##      0.22

``` r
# OPC algorithm with the ordinal Jonckheere–Terpstra test (significance level: 0.05)
OPCfit <- pc(suffStat = list(data = ordinal_data, stat_type = "jt"),
             alpha = 0.05,
             indepTest = catCItest,
             labels = colnames(ordinal_data))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
# Compare the patterns between them
comparePatterns(OPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     25.00      1.00      1.00    164.00     24.00      0.50      0.04      0.01 
    ##     FPR_P 
    ##      0.04

``` r
comparePatterns(OPCfit,trueDAG,hardP2P = FALSE) #soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     25.00      1.00      1.00    164.00     24.00      0.50      0.04      0.01 
    ##     FPR_P 
    ##      0.04

``` r
# GPC algorithm with the Gaussian test (significance level: 0.05)
GPCfit <- pc(suffStat = list(C = cor(ordinal_data), n = N),
             alpha = 0.05,
             indepTest = gaussCItest,
             labels = colnames(ordinal_data))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
# Compare the patterns between them
comparePatterns(GPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     29.00      9.00     21.00    152.00      8.00      0.30      0.36      0.13 
    ##     FPR_P 
    ##      0.84

``` r
comparePatterns(GPCfit,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     25.50     12.50     17.50    152.00      8.00      0.42      0.50      0.11 
    ##     FPR_P 
    ##      0.70

``` r
# hybrid method with the BDe score and the nominal PC output as the initial search space
BDE <- scoreparameters("bdecat",data.frame(ordinal_data),bdecatpar = list(chi = 0.5))
BDEfit <- iterativeMCMC(BDE)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
# Compare the patterns between them
comparePatterns(BDEfit$DAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     13.00     12.00      6.00    165.00      7.00      0.67      0.48      0.04 
    ##     FPR_P 
    ##      0.24

``` r
comparePatterns(BDEfit$DAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     10.00     15.00      3.00    165.00      7.00      0.83      0.60      0.02 
    ##     FPR_P 
    ##      0.12

``` r
# hybrid method with the BGe score and the Gaussian PC output as the initial search space
BGE <- scoreparameters("bge", ordinal_data, bgepar = list(am = 0.5))
BGEfit <- iterativeMCMC(BGE,scoreout = TRUE)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
# Compare the patterns between them
comparePatterns(BGEfit$DAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     16.00     20.00     13.00    154.00      3.00      0.61      0.80      0.08 
    ##     FPR_P 
    ##      0.52

``` r
comparePatterns(BGEfit$DAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     15.50     20.50     12.50    154.00      3.00      0.62      0.82      0.08 
    ##     FPR_P 
    ##      0.50

``` r
# the OSEM algorithm
OSEMfit <- ordinalStructEM(n, ordinal_data,
                           usrpar = list(penType = "other",
                                         L = 5,
                                         lambda = 2))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-25-1.png)

``` r
# Compare the patterns between them
comparePatterns(OSEMfit$DAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##      7.00     24.00      7.00    159.00      0.00      0.77      0.96      0.04 
    ##     FPR_P 
    ##      0.28

``` r
comparePatterns(OSEMfit$DAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##      6.50     24.50      6.50    159.00      0.00      0.79      0.98      0.04 
    ##     FPR_P 
    ##      0.26

(Tsagris et al., 2018)

``` r
d <- as.data.frame(ordinal_data)
d[] <- lapply(d[], as.ordered)
skel <- pc.skel(d, method = "comb.mm", alpha = 0.05)
MMDAG <- pc.or(skel)$G
MMDAG[MMDAG == 2] <- 1
MMDAG[MMDAG == 3] <- 0
```

![](Demo_files/figure-markdown_github/unnamed-chunk-28-1.png)

``` r
# Compare the patterns between them
comparePatterns(MMDAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     25.00      8.00     17.00    157.00      8.00      0.32      0.32      0.10 
    ##     FPR_P 
    ##      0.68

``` r
comparePatterns(MMDAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     22.00     11.00     14.00    157.00      8.00      0.44      0.44      0.08 
    ##     FPR_P 
    ##      0.56

(Cui et al., 2016)

Rank PC:

``` r
corr.rank <- sin(pi/2 * cor(ordinal_data, use = 'pairwise.complete.obs', method = 'kendall'))
RPCfit <- pc(suffStat = list(C = corr.rank, n = N), 
             indepTest = gaussCItest, labels = colnames(ordinal_data), alpha = 0.05, conservative = T)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-31-1.png)

``` r
# Compare the patterns between them
comparePatterns(RPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     23.00     15.00     17.00    152.00      6.00      0.47      0.60      0.10 
    ##     FPR_P 
    ##      0.68

``` r
comparePatterns(RPCfit,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     21.50     16.50     15.50    152.00      6.00      0.52      0.66      0.09 
    ##     FPR_P 
    ##      0.62

Copula PC:

``` r
# copula object
cop.obj <- inferCopulaModel(ordinal_data, nsamp = 1000, S0 = diag(n)/N, verb = F)
# correlation matrix samples
C_samples <- cop.obj$C.psamp[,, 501:1000]
# average correlation matrix
corr.cop <- apply(C_samples, c(1,2), mean)
# call the PC algorithm for causal discovery
CPCfit <- pc(suffStat = list(C = corr.cop, n = N), 
             indepTest = gaussCItest, labels = colnames(ordinal_data), alpha = 0.05, conservative = T)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-34-1.png)

``` r
# Compare the patterns between them
comparePatterns(CPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     19.00     16.00     12.00    155.00      7.00      0.57      0.64      0.07 
    ##     FPR_P 
    ##      0.48

``` r
comparePatterns(CPCfit,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     18.00     17.00     11.00    155.00      7.00      0.61      0.68      0.07 
    ##     FPR_P 
    ##      0.44
