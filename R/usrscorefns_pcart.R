require(rpcart)

### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam, usrpar){
  
  initparam$preLevels <- usrpar$preLevels
  initparam$pcart_alpha <- usrpar$pcart_alpha
  initparam$pcart_kappa <- usrpar$pcart_kappa
  initparam$df <- initparam$data
  initparam$data <- apply(as.matrix(initparam$data),2,function (x) as.numeric(x))
  return(initparam)
  
}

### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j,parentnodes,n,param) {
  
  res <- opt.pcart(param$data, parentnodes, j, param$preLevels, 
                   alpha = param$pcart_alpha, kappa = param$pcart_kappa)
  
  return(res$dataScore + res$structureScore)
}


# to.factor <- function(x) {
#   # Ensures that missing levels are not removed
#   if(class(x) == "factor") {
#     return(x)
#   } else {
#     return(factor(x))
#   }
# }