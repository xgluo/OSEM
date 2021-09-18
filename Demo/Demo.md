### 1. Load necessary packages and source code

``` r
library(BiDAG)
library(Rgraphviz)
library(pcalg)
library(graph)
library(bnlearn)
library(MXM)
library(sbgcop)
library(infotheo)
library(rpcart)
library(devtools)
source_url("https://raw.githubusercontent.com/cuiruifei/CausalMissingValues/master/R/inferCopulaModel.R")
```

``` r
setwd("../R")
# Major file containing the OSEM algorithm
source("ordinalScore.R") 
# Modify some of the existing functions in the BiDAG package to accommodate our user-defined functions
insertSource("spacefns.R",package = "BiDAG")
insertSource("usrscorefns.R",package = "BiDAG")
insertSource("initpar.R",package = "BiDAG")
insertSource("scoreagainstdag.R",package = "BiDAG")
```

### 2. Generate DAGs and ordinal data

``` r
set.seed(222)
# Generate a regular DAG with 20 nodes with 4 number of neighbors
n <- 20
trueDAG <- randDAG(n = n, d = 4, method = "er", wFUN = list(mywFUN))

# Generate a Gaussian dataset based on the DAG and standardize each dimension
# Sample size = 500
N <- 500
hidden_data <- rmvDAG2(N, trueDAG)
scaled_data <- t(t(hidden_data) - apply(hidden_data,2,mean))
truecov <- trueCov(trueDAG)
D <- diag(sqrt(diag(truecov)))
D.inv <- chol2inv(chol(D))
trueSigma <- D.inv %*% truecov %*% D.inv
scaled_data <- t(D.inv %*% t(scaled_data))

# Convert the Gaussian dataset into an ordinal dataset
ordinal_data <- convertToOrdinal(scaled_data, exp_levels = 4,concent_param = 2)
ordinal_data_df <- as.data.frame(ordinal_data)
ordinal_data_df[] <- lapply(ordinal_data_df[], as.ordered)
ordinal_levels <- apply(ordinal_data, 2, function(x) length(unique(x)))
```

### 3. Learn the structures using different methods

The following simulations can be repeated many times for different
configurations of parameters in order to obtain the ROC curves.

-   NPC

``` r
# NPC algorithm with the G^2 test (significance level: 0.05)
NPCfit <- pc(suffStat = list(dm = ordinal_data, 
                             nlev = apply(ordinal_data, 2, function (x) length(unique(x))),
                             adaptDF = FALSE),
             alpha = 0.05,
             indepTest = disCItest,
             labels = colnames(ordinal_data))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
# Compare the patterns between them
comparePatterns(NPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     37.50      2.50      4.50    150.00     33.00      0.36      0.06      0.03 
    ##     FPR_P 
    ##      0.12

``` r
comparePatterns(NPCfit,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     37.50      2.50      4.50    150.00     33.00      0.36      0.06      0.03 
    ##     FPR_P 
    ##      0.12

-   OPC (Musella, 2013)

``` r
# OPC algorithm with the ordinal Jonckheere–Terpstra test (significance level: 0.05) (bnlearn implementation)
OPCfit <- amat(pc.stable(ordinal_data_df, alpha = 0.05, test = "jt"))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# Compare the patterns between them
comparePatterns(OPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     29.50     14.50     18.50    146.00     11.00      0.44      0.37      0.12 
    ##     FPR_P 
    ##      0.47

``` r
comparePatterns(OPCfit,trueDAG,hardP2P = FALSE) #soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     29.50     14.50     18.50    146.00     11.00      0.44      0.37      0.12 
    ##     FPR_P 
    ##      0.47

-   GPC

``` r
# GPC algorithm with the Gaussian test (significance level: 0.05)
GPCfit <- pc(suffStat = list(C = cor(ordinal_data), n = N),
             alpha = 0.05,
             indepTest = gaussCItest,
             labels = colnames(ordinal_data))
# # GPC algorithm (significance level: 0.05) (bnlearn implementation)
# GPCfit <- amat(pc.stable(as.data.frame(ordinal_data), alpha = 0.05, test = "zf"))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-7-1.png)

``` r
# Compare the patterns between them
comparePatterns(GPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     30.00     13.00     21.00    147.00      9.00      0.38      0.33      0.14 
    ##     FPR_P 
    ##      0.54

``` r
comparePatterns(GPCfit,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     30.00     13.00     21.00    147.00      9.00      0.38      0.33      0.14 
    ##     FPR_P 
    ##      0.54

-   RPC (Harris and Drton, 2013; Cui et al., 2018)

``` r
corr.rank <- sin(pi/2 * cor(ordinal_data, use = 'pairwise.complete.obs', method = 'kendall'))
RPCfit <- pc(suffStat = list(C = corr.rank, n = N), 
             indepTest = gaussCItest, labels = colnames(ordinal_data), alpha = 0.05, conservative = T)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
# Compare the patterns between them
comparePatterns(RPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     17.00     24.00      7.00    149.00     10.00      0.77      0.62      0.05 
    ##     FPR_P 
    ##      0.18

``` r
comparePatterns(RPCfit,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     17.00     24.00      7.00    149.00     10.00      0.77      0.62      0.05 
    ##     FPR_P 
    ##      0.18

-   Copula PC (Cui et al., 2016, 2017, 2018)

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

![](Demo_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
# Compare the patterns between them
comparePatterns(CPCfit,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     28.00     16.00     13.00    146.00     15.00      0.55      0.41      0.09 
    ##     FPR_P 
    ##      0.33

``` r
comparePatterns(CPCfit,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     28.00     16.00     13.00    146.00     15.00      0.55      0.41      0.09 
    ##     FPR_P 
    ##      0.33

-   MMPC (Tsagris et al., 2018)

``` r
skel <- pc.skel(ordinal_data_df, method = "comb.mm", alpha = 0.05)
MMDAG <- pc.or(skel)$G
MMDAG[MMDAG == 2] <- 1
MMDAG[MMDAG == 3] <- 0
```

![](Demo_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
# Compare the patterns between them
comparePatterns(MMDAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     29.50     14.50     17.50    146.00     12.00      0.45      0.37      0.12 
    ##     FPR_P 
    ##      0.45

``` r
comparePatterns(MMDAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     29.50     14.50     17.50    146.00     12.00      0.45      0.37      0.12 
    ##     FPR_P 
    ##      0.45

-   BDe (Heckerman and Geiger, 1995)

``` r
# hybrid method with the BDe score and the nominal PC output as the initial search space
BDE <- scoreparameters("bdecat",data.frame(ordinal_data),bdecatpar = list(chi = 0.5))
BDEfit <- iterativeMCMC(BDE)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
# Compare the patterns between them
comparePatterns(BDEfit$DAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     26.00     13.00      8.00    151.00     18.00      0.62      0.33      0.05 
    ##     FPR_P 
    ##      0.21

``` r
comparePatterns(BDEfit$DAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     26.00     13.00      8.00    151.00     18.00      0.62      0.33      0.05 
    ##     FPR_P 
    ##      0.21

-   BGe (Heckerman and Geiger, 1995)

``` r
# hybrid method with the BGe score and the Gaussian PC output as the initial search space
BGE <- scoreparameters("bge", ordinal_data, bgepar = list(am = 0.5))
BGEfit <- iterativeMCMC(BGE,scoreout = TRUE)
```

![](Demo_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
# Compare the patterns between them
comparePatterns(BGEfit$DAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     13.00     35.00     10.00    142.00      3.00      0.78      0.90      0.07 
    ##     FPR_P 
    ##      0.26

``` r
comparePatterns(BGEfit$DAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     13.00     35.00     10.00    142.00      3.00      0.78      0.90      0.07 
    ##     FPR_P 
    ##      0.26

-   OSEM

``` r
# the OSEM algorithm
OSEMfit <- ordinalStructEM(n, ordinal_data,
                           usrpar = list(penType = "other",
                                         L = 5,
                                         lambda = 3,
                                         preLevels = NULL))
```

![](Demo_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
# Compare the patterns between them
comparePatterns(OSEMfit$DAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##      3.00     38.00      2.00    149.00      1.00      0.95      0.97      0.01 
    ##     FPR_P 
    ##      0.05

``` r
comparePatterns(OSEMfit$DAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##      3.00     38.00      2.00    149.00      1.00      0.95      0.97      0.01 
    ##     FPR_P 
    ##      0.05

-   PCART (Talvitie et al., 2019)

``` r
setwd("../R")
insertSource("usrscorefns_pcart.R",package = "BiDAG")
```

``` r
pcartparam <- scoreparameters("usr", ordinal_data_df, 
                              usrpar = list(pcart_alpha = 1.5,
                                            pcart_kappa = 0.25,
                                            pctesttype = "bde", 
                                            preLevels = ordinal_levels,
                                            response_type = "CAT"))
# need to set limit to the parent set size due to computational limit
pcartfit <- iterativeMCMC(pcartparam, alpha = 0, plus1it = 10, softlimit = 3, hardlimit = 3)
```

    ## maximum parent set size is 0 
    ## core space defined, score table are being computed 
    ## score tables completed, iterative MCMC is running 
    ## search space expansion 2 
    ## search space expansion 3 
    ## search space expansion 4 
    ## search space expansion 5

![](Demo_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
# Compare the patterns between them
comparePatterns(pcartfit$DAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     26.00     17.00     16.00    147.00     10.00      0.52      0.44      0.11 
    ##     FPR_P 
    ##      0.41

``` r
comparePatterns(pcartfit$DAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     26.00     17.00     16.00    147.00     10.00      0.52      0.44      0.11 
    ##     FPR_P 
    ##      0.41

-   OPCART (Talvitie et al., 2019)

``` r
opcartparam <- scoreparameters("usr", ordinal_data_df, 
                              usrpar = list(pcart_alpha = 1.5,
                                            pcart_kappa = 0.25,
                                            pctesttype = "bde", 
                                            preLevels = ordinal_levels,
                                            response_type = "ORD"))
# need to set limit to the parent set size due to computational limit
opcartfit <- iterativeMCMC(opcartparam, alpha = 0, plus1it = 10, softlimit = 3, hardlimit = 3)
```

    ## maximum parent set size is 0 
    ## core space defined, score table are being computed 
    ## score tables completed, iterative MCMC is running 
    ## search space expansion 2 
    ## search space expansion 3 
    ## search space expansion 4

![](Demo_files/figure-markdown_github/unnamed-chunk-24-1.png)

``` r
# Compare the patterns between them
comparePatterns(opcartfit$DAG,trueDAG) # hard version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     42.00      5.00      9.00    143.00     33.00      0.36      0.13      0.06 
    ##     FPR_P 
    ##      0.23

``` r
comparePatterns(opcartfit$DAG,trueDAG,hardP2P = FALSE) # soft version
```

    ##       SHD        TP        FP        TN        FN Precision       TPR     FPR_N 
    ##     42.00      5.00      9.00    143.00     33.00      0.36      0.13      0.06 
    ##     FPR_P 
    ##      0.23
