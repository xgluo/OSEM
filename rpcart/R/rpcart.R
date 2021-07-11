
##' @export
opt.pcart <- function(data, predictors, response, preLevels, response_type = "ORD", 
                      alpha = 0.5, kappa = 0.25, ess = 1.0, use_structure_score=TRUE) {
  
  rdata <- as.matrix(data[,c(predictors,response)])
  res <- cli_pcart(rdata,
                   length(predictors), nrow(data), preLevels[predictors],
                   response_type, preLevels[response],
                   alpha, kappa, ess, use_structure_score)
  
  return(res)
}
