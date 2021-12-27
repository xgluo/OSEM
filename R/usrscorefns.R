### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam,
                               usrpar = list(penType = c("AIC","BIC","other"),
                                             L = 5,
                                             lambda = 2,
                                             pctesttype = "usr")){
  
  return(ordinalScoreParam(initparam, usrpar))
  
}

### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j,parentnodes,n,param) {
  
  return(ordinalCoreScore(j,parentnodes,n,param))
}

### This function defines the CI tests for the starting skeleton
usrdefinestartspace <- function(alpha,param,cpdag,n){

  cormat<-param$Sigma_hat
  N <- param$N
  if(cpdag){
    pc.skel<-pcalg::pc(suffStat = list(C = cormat, n = N),
                       indepTest = gaussCItest,
                       alpha=alpha,labels=colnames(param$data),skel.method="stable",verbose = FALSE)
  } else {
    pc.skel<-pcalg::skeleton(suffStat = list(C = cormat, n = N),
                             indepTest = gaussCItest,
                             alpha=alpha,labels=colnames(param$data),method="stable",verbose = FALSE)
  }
  
  return(pc.skel)
}
