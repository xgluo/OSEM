newspaceskel<-function(n,startspace,currspace,softlimit,hardlimit,posterior,blacklist,
                       MCMCtrace=NULL,mergetype="skeleton") {
  
  switch(mergetype,
         "dag" = { 
           mdag<-modelpcore(MCMCtrace,p=posterior,pdag=FALSE)
           newadj<-1*(!blacklist&(startspace|mdag))
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-(1*(!blacklist&mdag))[,toomanyneib]}
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-currspace[,toomanyneib]}
         },
         "cpdag" = { 
           mcp<-modelpcore(n,MCMCtrace,p=posterior,pdag=TRUE)
           newadj<-1*(!blacklist&(startspace|mcp))
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             mdag<-modelpcore(MCMCtrace,p=posterior,pdag=FALSE)
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|mdag)))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-mdag[,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
         },
         "skeleton" = { 
           mskel<-1*(modelpcore(MCMCtrace,p=posterior,pdag=FALSE)|t(modelpcore(MCMCtrace,p=posterior,pdag=FALSE)))
           newadj<-1*(!blacklist&(startspace|mskel))
           toomanyneib<-which(apply(newadj,2,sum)>4)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|modelpcore(MCMCtrace,p=posterior,pdag=TRUE))))[,toomanyneib]
           }
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             mdag<-modelpcore(MCMCtrace,p=posterior,pdag=FALSE)
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|mdag)))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-mdag[,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
         }
  )
  return(newadj)
}

newspacemap<-function(n,startspace,currspace,softlimit,hardlimit,blacklist, 
                      maxdag=NULL,mergetype="skeleton",accum) {
  
  if(!is.matrix(maxdag)) maxdag<-as.matrix(maxdag)
  switch(mergetype,
         "dag" = { 
           maxdag<-maxdag
           newadj<-1*(!blacklist&(startspace|maxdag))
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-(1*(!blacklist&maxdag))[,toomanyneib]}
           toomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(toomanyneib)>0){newadj[,toomanyneib]<-currspace[,toomanyneib]}
           if(accum) {
             newnewadj<-1*(newadj|currspace)
             toomanyneib<-which(apply(newnewadj,2,sum)>hardlimit)
             if(length(toomanyneib)>0){newnewadj[,toomanyneib]<-newadj[,toomanyneib]}
             newadj<-newnewadj
           }
         },
         "cpdag" = { 
           maxcp<-dagadj2cpadj(maxdag)
           newadj<-1*(!blacklist&(startspace|maxcp))
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|maxdag)))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-maxdag[,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
           if(accum) {
             newnewadj<-1*(newadj|currspace)
             toomanyneib<-which(apply(newnewadj,2,sum)>hardlimit)
             if(length(toomanyneib)>0){newnewadj[,toomanyneib]<-newadj[,toomanyneib]}
             newadj<-newnewadj
           }
         },
         "skeleton" = { 
           maxskel<-1*(maxdag|transp(maxdag))
           newadj<-1*(!blacklist&(startspace|maxskel))
           toomanyneib<-which(apply(newadj,2,sum)>7)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|dagadj2cpadj(maxdag))))[,toomanyneib]
           }
           toomanyneib<-which(apply(newadj,2,sum)>softlimit)
           if(length(toomanyneib)>0) {
             newadj[,toomanyneib]<-(1*(!blacklist&(startspace|maxdag)))[,toomanyneib]
           }
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-maxdag[,tootoomanyneib]}
           tootoomanyneib<-which(apply(newadj,2,sum)>hardlimit)
           if(length(tootoomanyneib)>0) {newadj[,tootoomanyneib]<-currspace[,tootoomanyneib]}
           if(accum) {
             newnewadj<-1*(newadj|currspace)
             toomanyneib<-which(apply(newnewadj,2,sum)>hardlimit)
             if(length(toomanyneib)>0){newnewadj[,toomanyneib]<-newadj[,toomanyneib]}
             newadj<-newnewadj
           }
         }
  )
  return(newadj)
}


definestartspace<-function(alpha,param,cpdag=FALSE,algo="pc",alphainit=NULL) {
  if(is.null(alphainit)) {alphainit<-alpha}
  
  local_type <- param$type
  if(local_type=="usr") {
    if(param$pctesttype%in%c("bde","bge","bdecat")) {
      local_type<-param$pctesttype
    } 
  }
  
  if(param$DBN){
    if(param$stationary) {
      othersliceskel <- definestartspace(alpha,param$otherslices,cpdag=FALSE,algo="pc")
      firstsliceskel <- definestartspace(alphainit,param$firstslice,cpdag=FALSE,algo="pc")
      startspace <- othersliceskel
      startspace[param$intstr$rows,param$intstr$cols] <- 1*(startspace[param$intstr$rows,param$intstr$cols] | firstsliceskel[param$intstr$rows,param$intstr$cols])
      #diag(startspace[param$trans$rows,param$trans$cols])<-1
    } else {
      skels<-list()
      skels[[1]]<-definestartspace(alphainit,param$paramsets[[1]],cpdag=FALSE,algo="pc")
      startspace<-skels[[1]]
      for(i in 2:(length(param$paramsets)-1)) {
        skels[[i]]<-definestartspace(alpha,param$paramsets[[i]],cpdag=FALSE,algo="pc")
        startspace<-1*(skels[[i]]|startspace)
      }
      firstsliceskel <- definestartspace(alphainit,param$paramsets[[length(param$paramsets)]],cpdag=FALSE,algo="pc")
      startspace[param$intstr$rows,param$intstr$cols] <- 1*(startspace[param$intstr$rows,param$intstr$cols] | firstsliceskel[param$intstr$rows,param$intstr$cols])
      #diag(startspace[param$trans$rows,param$trans$cols])<-1
    }
  } else { # otherwise use old versions
    
    if(local_type=="bde") {
      if(cpdag){
        pc.skel<-pc(suffStat = list(d1=param$d1,d0=param$d0,data=param$data),
                    indepTest = weightedbinCItest, alpha = alpha, labels = colnames(param$data),
                    verbose = FALSE)
        
      } else {
        pc.skel<-pcalg::skeleton(suffStat = list(d1=param$d1,d0=param$d0,data=param$data),
                                 indepTest = weightedbinCItest, alpha = alpha, labels = colnames(param$data),
                                 verbose = FALSE)
      }
    } else if(local_type=="bdecat") {
      if(cpdag){
        pc.skel<-pc(suffStat = param,
                    indepTest = weightedcatCItest, alpha = alpha, labels = colnames(param$data),
                    verbose = FALSE)
        
      } else {
        pc.skel<-pcalg::skeleton(suffStat = param,
                                 indepTest = weightedcatCItest, alpha = alpha, labels = colnames(param$data),
                                 verbose = FALSE)
      }
    } else if(local_type=="bge") {
      if(is.null(param$weightvector)) {
        cormat<-cor(param$data)
        N<-nrow(param$data)
      } else { N<-sum(param$weightvector)
      cormat<-cov.wt(param$data,wt=param$weightvector,cor=TRUE)$cor}
      if(cpdag){
        pc.skel<-pcalg::pc(suffStat = list(C = cormat, n = N),
                           indepTest = gaussCItest,
                           alpha=alpha,labels=colnames(param$data),skel.method="stable",verbose = FALSE)
      } else {
        pc.skel<-pcalg::skeleton(suffStat = list(C = cormat, n = N),
                                 indepTest = gaussCItest,
                                 alpha=alpha,labels=colnames(param$data),method="stable",verbose = FALSE)
      }
    } else if (local_type=="usr") {
      
      pc.skel <- usrdefinestartspace(alpha,param,cpdag,n)
    }
    
    g<-pc.skel@graph
    startspace<-1*(graph2m(g))
  }
  return(startspace)
}

