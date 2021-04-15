require(rpcart)

### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam, usrpar){
  
  return(initparam)
  
}

### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j,parentnodes,n,param) {
  
  res <- opt.pcart.cat(param$data, parentnodes, j)
  
  return(res$score)
}


