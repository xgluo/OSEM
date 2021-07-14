#'Estimating posterior probabilities of single edges
#'
#'This function estimates the posterior probabilities of edges by averaging over a sample of DAGs
#'obtained via an MCMC scheme.
#'
#'@param MCMCchain an object of class \code{partitionMCMC} or \code{orderMCMC}, representing the output of structure sampling function \code{\link{partitionMCMC}} or \code{\link{orderMCMC}} (the latter when parameter \code{chainout}=TRUE; 
#'@param pdag logical, if TRUE (FALSE by default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin number between \code{0} and \code{1}, indicates the percentage of the samples which will be discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@param endstep number between \code{0} and \code{1}; 1 by default 
#'@return a square matrix with dimensions equal to the number of variables; each entry \code{[i,j]} is an estimate of the posterior probability of the edge from node \code{i} to node \code{j}
#'@examples
#'Bostonscore<-scoreparameters("bge", Boston)
#'\dontrun{
#'samplefit<-orderMCMC(Bostonscore, iterations=25000,chainout=TRUE)
#'edgesposterior<-edgep(samplefit, pdag=TRUE, burnin=0.2)
#'}
#'@author Polina Suter
#'@export
edgep<-function(MCMCchain,pdag=FALSE,burnin=0.2,endstep=1) {
  
  DBN<-MCMCchain$info$DBN
  MCMCinfo<-MCMCchain$info
  MCMCchain<-MCMCchain$traceadd$incidence
    
  varlabels<-colnames(MCMCchain[[1]])
  if(endstep==1) {
  endstep<-length(MCMCchain)
  } else {
    endstep<-ceiling(length(MCMCchain)*endstep)
  }
  startstep<-max(as.integer(burnin*endstep),1)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    incidence<-Reduce('+', cpdags)/(endstep-startstep+1)
  } else {
    incidence<-Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1)
  }
  colnames(incidence)<-varlabels
  rownames(incidence)<-varlabels
  if(DBN) {
    incidence<-DBNcut(incidence,n.dynamic=MCMCinfo$nsmall,n.static=MCMCinfo$bgn)
    incidence.init<-DBNinit(incidence,n.dynamic=MCMCinfo$nsmall,n.static=MCMCinfo$bgn)
    incidence[1:(MCMCinfo$nsmall+MCMCinfo$bgn),1:(MCMCinfo$nsmall+MCMCinfo$bgn)]<-incidence.init
  }
  return(incidence)
}


#'Estimating a graph corresponding to a posterior probability threshold
#'
#'This function constructs a directed graph (not necessarily acyclic) including all edges with a posterior probability above a certain threshold.  The posterior probability is evaluated as the Monte Carlo estimate from a sample of DAGs obtained via an MCMC scheme.
#'
#'@param MCMCchain object of class \code{partitionMCMC} or \code{orderMCMC}, representing the output of structure sampling function \code{\link{partitionMCMC}} or \code{\link{orderMCMC}} (the latter when parameter \code{chainout}=TRUE; 
#'@param p threshold such that only edges with a higher posterior probability will be retained in the directed graph summarising the sample of DAGs
#'@param pdag logical, if TRUE (FALSE by default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin number between \code{0} and \code{1}, indicates the percentage of the samples which will be  the discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@return a square matrix with dimensions equal to the number of variables representing the adjacency matrix of the directed graph summarising the sample of DAGs
#'@examples
#'Bostonscore<-scoreparameters("bge", Boston)
#'\dontrun{
#'orderfit<-orderMCMC(Bostonscore, MAP=FALSE, iterations=25000, chainout=TRUE)
#'hdag<-modelp(orderfit, p=0.9)
#'}
#'@author Polina Suter
#'@export
modelp<-function(MCMCchain, p, pdag=FALSE, burnin=0.2) {

  DBN<-MCMCchain$info$DBN
  MCMCinfo<-MCMCchain$info
  MCMCchain<-MCMCchain$traceadd$incidence

  varlabels<-colnames(MCMCchain[[1]])
  n<-nrow(MCMCchain[[1]])
  incidence<-matrix(rep(0, n*n), nrow=n, ncol=n)
  endstep<-length(MCMCchain)
  startstep<-max(as.integer(burnin*endstep),1)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    incidence[which(Reduce('+', cpdags)/(endstep-startstep+1)>p)]<-1
  } else {
    incidence[which(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1)>p)]<-1
  }
  colnames(incidence)<-varlabels
  rownames(incidence)<-varlabels
  if(DBN) {
      incidence<-DBNcut(incidence,n.dynamic=MCMCinfo$nsmall,n.static=MCMCinfo$bgn)
      incidence.init<-DBNinit(incidence,n.dynamic=MCMCinfo$nsmall,n.static=MCMCinfo$bgn)
      incidence[1:(MCMCinfo$nsmall+MCMCinfo$bgn),1:(MCMCinfo$nsmall+MCMCinfo$bgn)]<-incidence.init
  }
  return(incidence)
}


#'Performance assessment of iterative MCMC scheme against a known Bayesian network
#'
#'This function compute 8 different metrics of structure fit of an object of class \code{iterativeMCMC} to the ground truth DAG (or CPDAG). Object of class
#'\code{iterativeMCMC} stores MAP graph at from each search space expansion step. This function computes structure fit of
#'each of the stored graphs to the ground truth one. Computed metrics include: TP, FP, TPR, FPR, FPRn, FDR, SHD. See metrics description in
#'see also \code{\link{compareDAGs}}.
#'
#'@param MCMCmult an object which of class \code{iterativeMCMC}, see also \code{\link{iterativeMCMC}})
#'@param truedag ground truth DAG which generated the data used in the search procedure; represented by an object of class \code{\link[graph]{graphNEL}} or an adjacency matrix
#@param sample logical (FALSE by default), indicates if \code{MCMCmult} contains sample or maximum score DAGs
#'@param cpdag logical, if TRUE (FALSE by default) all DAGs are first converted to their respective equivalence classes (CPDAG) 
#'@param p threshold such that only edges with a higher posterior probability will be retained in the directed graph summarising the sample of DAGs at each iteration from \code{MCMCmult} if parameter \code{sample} set to TRUE
#'@param trans logical, for DBNs indicates if model comparions are performed for transition structure; when \code{trans} equals FALSE the comparison is performed for initial structures of estimated models and the ground truth DBN; for usual BNs the parameter is disregarded
#'@return an object if class \code{itersim}, a matrix with the number of rows equal to the number of expansion iterations in \code{iterativeMCMC}, and 8 columns reporting for 
#'the maximally scoring DAG uncovered at each iteration: the number of true positive edges ('TP'), the number of false positive edges ('FP'), 
#'the true positive rate ('TPR'), the structural Hamming distance ('SHD'), false positive rate ('FPR'),
#'false discovery rate ('FDR') and the score of the DAG (`score'). 
#' @examples
#' gsim.score<-scoreparameters("bge", gsim)
#' \dontrun{
#' MAPestimate<-iterativeMCMC(gsim.score)
#' itercomp(MAPestimate, gsimmat)
#' }
#'@author Polina Suter
#'@export
itercomp<-function(MCMCmult, truedag, cpdag=TRUE, p=0.5,trans=TRUE) {

  if(MCMCmult$info$DBN) { #we need to extract either transition or initial structure
     
       if(trans) {
        if(!is.matrix(truedag)) truedag<-graph2m(truedag)
        truedag<-m2graph(DBNcut(truedag,n.dynamic=MCMCmult$info$nsmall,n.static=MCMCmult$info$bgn))
        trueskeleton<-graph2skeleton(truedag)
      } else {
        if(!is.matrix(truedag)) truedag<-graph2m(truedag)
        truedag<-m2graph(DBNinit(truedag,n.dynamic=MCMCmult$info$nsmall,n.static=MCMCmult$info$bgn))
        trueskeleton<-graph2skeleton(truedag)
      }
      
      if(MCMCmult$info$split) {
        if(trans) {
          MCMCmult$maxtrace<-MCMCmult$trans$maxtrace
          MCMCmult$trans$maxtrace<-NULL
        } else {
          MCMCmult$maxtrace<-MCMCmult$init$maxtrace
          MCMCmult$init$maxtrace<-NULL
        }
      } else {
          newtrace<-lapply(MCMCmult$maxtrace,function(x)x$DAG)
          if(trans) {
            newtrace<-lapply(newtrace,DBNcut,n.dynamic=MCMCmult$info$nsmall,n.static=MCMCmult$info$bgn)
            for(i in 1:length(newtrace)) {
              MCMCmult$maxtrace[[i]]$DAG<-newtrace[[i]]
            }
          } else {
            newtrace<-lapply(newtrace,DBNinit,n.dynamic=MCMCmult$info$nsmall,n.static=MCMCmult$info$bgn)
            for(i in 1:length(newtrace)) {
              MCMCmult$maxtrace[[i]]$DAG<-newtrace[[i]]
            }
            }
      }
    } 
  
  mapdags<-lapply(MCMCmult$maxtrace, function(x)x$DAG)
  score<-unlist(lapply(MCMCmult$maxtrace, function(x)x$score))
  res<-Reduce(rbind, lapply(mapdags, compareDAGs, truedag, cpdag=cpdag))
  res<-cbind(res,score)
  rownames(res)<-c(1:nrow(res))
  attr(res,"class")<-"itercomp"
  return(res) 
}


#'Performance assessment of sampling algorithms against a known Bayesian network
#'
#'This function compute 8 different metrics of structure fit of an object of classes \code{orderMCMC} and \code{partitionMCMC} to the ground truth DAG (or CPDAG). First posterior probabilities
#'of single edges are calculated based on a sample stores in the object of class \code{orderMCMC} or \code{partitionMCMC}. This function computes structure fit of
#'each of the consensus graphs to the ground truth one based on a defined range of posterior thresholds. Computed metrics include: TP, FP, TPR, FPR, FPRn, FDR, SHD. See metrics description in
#'see also \code{\link{compareDAGs}}.
#'
#'@param MCMCchain an object of class \code{partitionMCMC} or \code{orderMCMC}, representing the output of structure sampling function \code{\link{partitionMCMC}} or \code{\link{orderMCMC}} (the latter when parameter \code{chainout}=TRUE; 
#'@param truedag ground truth DAG which generated the data used in the search procedure; represented by an object of class  \code{\link[graph]{graphNEL}}
#'@param p a vector of numeric values between 0 and 1, defining posterior probabilities according to which the edges of assessed structures are drawn, please note very low barriers can lead to very dense structures; by default 
#'\eqn{p=c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)}
#'@param pdag logical, if TRUE (default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging
#'@param burnin number between \code{0} and \code{1}, indicates the percentage of the samples which will be  the discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default
#'@param trans logical, for DBNs indicates if model comparions are performed for transition structure; when \code{trans} equals FALSE the comparison is performed for initial structures of estimated models and the ground truth DBN; for usual BNs the parameter is disregarded
#'@return an object if class \code{samplesim}, a matrix with the number of rows equal to the number of elements in 'p', and 8 columns reporting for 
#'the consensus graphss (corresponfing to each of the values in 'p') the number of true positive edges ('TP'), the number of false positive edges ('FP'), the number of false negative edges ('FN'),
#'the true positive rate ('TPR'), the structural Hamming distance ('SHD'), false positive rate ('FPR'),
#'false discovery rate ('FDR') and false positive rate normalized by TP+FN ('FPRn').
#' @examples
#' gsim.score<-scoreparameters("bge", gsim)
#' \dontrun{
#' mapest<-iterativeMCMC(gsim.score)
#' ordersample<-orderMCMC(gsim.score, MAP=FALSE, startspace=mapest$endspace)
#' samplecomp(ordersample, gsimmat)
#' }
#'@author Polina Suter
#'@export
samplecomp<-function(MCMCchain, truedag, p=c(0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2),
                                         pdag=TRUE, burnin=0.2, trans=TRUE) {
  if(is.matrix(truedag)) truedag<-m2graph(truedag)
  MCMCmatlist<-MCMCchain$traceadd$incidence
  n<-nrow(MCMCmatlist[[1]])
  truecp<-pcalg::dag2cpdag(truedag)
  if(MCMCchain$info$DBN) {
    pdag<-FALSE
    if(trans==TRUE) {
        cat("comparison is performed for transition structures \n")
        trueadj<-DBNcut(graph2m(truedag),MCMCchain$info$nsmall,MCMCchain$info$bgn)
        truedag<-m2graph(trueadj)
        truecpadj<-graph2m(truecp)
        truecpadj<-DBNcut(truecpadj,MCMCchain$info$nsmall,MCMCchain$info$bgn)
        truecp<-m2graph(truecpadj) 
      } else {
          cat("comparison is performed for initial structures \n")
          trueadj<-DBNinit(graph2m(truedag),MCMCchain$info$nsmall,MCMCchain$info$bgn)
          truedag<-m2graph(trueadj)
          truecpadj<-DBNinit(graph2m(truecp),MCMCchain$info$nsmall,MCMCchain$info$bgn)
          truecp<-m2graph(truecpadj)
          n<-MCMCchain$info$nsmall+MCMCchain$info$bgn
      }
  }

  endstep<-length(MCMCmatlist)
  startstep<-max(as.integer(burnin*endstep),1)
  
  if(pdag) {
  dags<-lapply(MCMCmatlist[startstep:endstep],dagadj2cpadj) #first convert every DAG in the sample to equivalence class
  } else {dags<-MCMCmatlist[startstep:endstep]}
  
  if(MCMCchain$info$DBN){
    if(trans) {
    dags<-lapply(dags,DBNcut,n.dynamic=MCMCchain$info$nsmall,n.static=MCMCchain$info$bgn)
    } else {
      dags<-lapply(dags,DBNinit,n.dynamic=MCMCchain$info$nsmall,n.static=MCMCchain$info$bgn)
    }
  }
  postprobmat<-Reduce('+', dags)/(endstep-startstep+1)
 
  if(length(p)==1) {
    mlist<-matrix(0, nrow=n,ncol=n)
    mlist[which(postprobmat>p)]<-1
    res<-compareDAGs(mlist,truedag, cpdag=pdag) 
    res<-c(res,p)
    names(res)[length(res)]<-"p"
  } else {
  mlist<-list()
  i<-1
  for (py in 1:length(p)) {
    mlist[[i]]<-matrix(0, nrow=n,ncol=n)
    mlist[[i]][which(postprobmat>p[py])]<-1
    i<-i+1
  }
  
  res<-lapply(mlist,compareDAGs,truedag, cpdag=pdag)
  res<-Reduce(rbind,res)
  res<-cbind(res,p)
  rownames(res)<-c(1:nrow(res))
  }
  attr(res,"class")<-"samplecomp"
  return(res)
}

modelpcore<-function(MCMCchain, p, pdag=FALSE, burnin=0.2, DBN=FALSE, nsmall=0, n.dynamic=0, n.static=0) {
  
  varlabels<-colnames(MCMCchain[[1]])
  n<-nrow(MCMCchain[[1]])
  incidence<-matrix(rep(0, n*n), nrow=n, ncol=n)
  endstep<-length(MCMCchain)
  startstep<-max(as.integer(burnin*endstep),1)
  if (pdag) {
    cpdags<-lapply(MCMCchain[startstep:endstep],dagadj2cpadj)
    incidence[which(Reduce('+', cpdags)/(endstep-startstep+1)>p)]<-1
  } else {
    incidence[which(Reduce('+', MCMCchain[startstep:endstep])/(endstep-startstep+1)>p)]<-1
  }
  colnames(incidence)<-varlabels
  rownames(incidence)<-varlabels
  if(DBN) {
    incidence<-DBNcut(incidence,n.dynamic=n.dynamic,n.static=n.static)
    incidence.init<-DBNinit(incidence,n.dynamic=n.dynamic,n.static=n.static)
    incidence[1:(n.dynamic+n.static),1:(n.dynamic+n.static)]<-incidence.init
  }
  return(incidence)
}
