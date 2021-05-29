require(rpcart)

### This function returns the objects needed to evaluate the user defined score
usrscoreparameters <- function(initparam, usrpar = list(preLevels = NULL,
                                                        pcart_alpha = 0.5,
                                                        pcart_kappa = 0.5,
                                                        response_type = "CAT",
                                                        pctesttype = "bde")){
  
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
    initparam$preLevels <- apply(initparam$data, 2, function (x) length(unique(x)))
  } else {
    initparam$preLevels <- usrpar$preLevels
  }

  initparam$pcart_alpha <- usrpar$pcart_alpha
  initparam$pcart_kappa <- usrpar$pcart_kappa
  initparam$response_type <- usrpar$response_type
  initparam$df <- initparam$data
  initparam$data <- apply(as.matrix(initparam$data),2,function (x) as.numeric(x))
  return(initparam)
  
}

### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j,parentnodes,n,param) {
  
  res <- opt.pcart(param$data, parentnodes, j, param$preLevels, 
                   alpha = param$pcart_alpha, kappa = param$pcart_kappa,
                   response_type = param$response_type)
  
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