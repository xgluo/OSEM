library(parallel)

cl <- makeCluster(10)

clusterEvalQ(cl,{
  library(pcalg)
  library(BiDAG)
  library(gtools)
  library(abind)
  load("OCDRogers.RData")
  source("ordinalScore.R") 
  insertSource("spacefns.R",package = "BiDAG")
  insertSource("usrscorefns.R",package = "BiDAG")
  insertSource("initpar.R",package = "BiDAG")
})

clusterEvalQ(cl, {
  sim_once <- function() {
    n <- ncol(datRogers)
    N <- nrow(datRogers)
    boot.sample <- datRogers[sample(N,N,replace = TRUE),]
    BDE <- scoreparameters(n,"bdecat",boot.sample,bdecatpar = list(chi = 1.5))
    BDEfit <- iterativeMCMC(BDE)
    BGE <- scoreparameters(n, "bge", boot.sample,bgepar = list(am = 0.05))
    BGEfit <- iterativeMCMC(BGE)
    OSEMfit <- ordinalStructEM(n, datRogers,
                               usrpar = list(penType = "other",
                                             L = 5,
                                             lambda = 6))
    res <- abind(BDEfit$max$DAG,BGEfit$max$DAG,OSEMfit$max$DAG, along = 3)
    dimnames(res)[[3]] <- c("BDe","BGe","OSEM")
    return(res)
  }
})

results <- parLapply(cl,1:500,function(i) sim_once())
stopCluster(cl)
res <- simplify2array(results)
save(res, file = "OCDRogers_boot.RData")
