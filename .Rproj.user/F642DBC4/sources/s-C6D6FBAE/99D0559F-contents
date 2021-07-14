#'Deriving an adjacency matrix of a graph
#'
#'This function derives the adjacency matrix corresponding to a graph object
#'
#'@param g graph, object of class \code{\link[graph]{graphNEL}} (package `graph')
#'@return a square matrix whose dimensions are the number of nodes in the graph g, where element
#' \code{[i,j]} equals \code{1} if there is a directed edge from node \code{i} to node \code{j} in the graph \code{g},
#'  and \code{0} otherwise
#' @examples 
#' Asiagraph<-m2graph(Asiamat)
#' Asia.adj<-graph2m(Asiagraph)
#'@export
graph2m<-function(g) {
  l<-length(g@edgeL)
  adj<-matrix(rep(0,l*l),nrow=l,ncol=l)
  for (i in 1:l) {
    adj[g@edgeL[[i]]$edges,i]<-1}
  rownames(adj)<-g@nodes
  colnames(adj)<-g@nodes
  return(t(adj))
}

#'Deriving connected subgraph 
#'
#'This function derives an adjacency matrix of a subgraph whose nodes are connected to at least one other node in a graph
#'
#'@param adj square adjacency matrix with elements in \code{\{0,1\}}, representing a graph
#'@return adjacency matrix of a subgraph of graph represented by 'adj' whose nodes have at least one connection 
#'@examples 
#'dim(gsimmat) #full graph contains 100 nodes
#'gconn<-connectedSubGraph(gsimmat) #removing disconnected nodes
#'dim(gconn) #connected subgraph contains 93 nodes
#'@export
connectedSubGraph<-function(adj) {
  n<-ncol(adj)
  col2del<-vector()
  for(i in 1:n) {
    if(sum(adj[i,],adj[,i])==0) {
      col2del<-c(col2del,i)
    }
  }
  if(length(col2del)>0) {
    adj<-adj[,-col2del]
    adj<-adj[-col2del,]
  }  
  return(adj)
}

#'Deriving subgraph 
#'
#'This function derives an adjacency matrix of a subgraph based on the adjacency matrix of a full graph and a list of nodes
#'
#'@param adj square adjacency matrix with elements in \code{\{0,1\}}, representing a graph
#'@param nodes vector of node names of the subgraph; should be a subset of column names of 'adj'
#'@return adjacency matrix of a subgraph which includes all 'nodes'
#'@examples 
#'getSubGraph(Asiamat,c("E","B","D","X"))
#'@export
getSubGraph<-function(adj,nodes) {
  if(!all(nodes%in%colnames(adj))) {
    stop("some 'nodes' can't be found in column names of adj!")
  }
  adj<-adj[nodes,nodes]
  return(adj)
}


#'Deriving a graph from an adjacancy matrix
#'
#'This function derives a graph object corresponding to an adjacency matrix
#'
#'@param adj square adjacency matrix with elements in \code{\{0,1\}}, representing a graph
#'@param nodes (optional) labels of the nodes, \code{c(1:n)} are used by default
#'@return object of class \code{\link[graph]{graphNEL}} (package `graph'); if element \code{adj[i,j]} equals \code{1}, then there is a directed edge from node \code{i} to node \code{j} in the graph, and no edge otherwise
#'@examples 
#'m2graph(Asiamat)
#'@export
m2graph<-function(adj,nodes=NULL) {
  l<-ncol(adj)
  if (is.null(nodes)) {
    if (!all(is.character(colnames(adj)))) {
    V <- c(1:l)
    edL <- vector("list", length=l)
    names(edL) <- sapply(V,toString)
    } else {
      edL <- vector("list", length=l)
      names(edL) <- colnames(adj)
      V<-colnames(adj)
    }
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))}
  else {
    V <- nodes
    edL <- vector("list", length=l)
    names(edL) <- V
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))
  }
  gR <- new("graphNEL", nodes=sapply(V,toString), edgeL=edL,edgemode="directed")
  
  return(gR)
  
}

#'Comparing two graphs
#'
#'This function compares one (estimated) graph to another graph (true graph), returning a vector of 8 values: 
#'\itemize{
#' \item the number of true positive edges ('TP') is the number of edges in the skeleton of 'egraph' which are also present in the skeleton of 'truegraph'
#' \item the number of false positive edges ('FP') is the number of edges in the skeleton of 'egraph' which are absent in the skeleton of 'truegraph'
#' \item the number of fralse negative edges ('FN') is the number of edges in the skeleton of 'truegraph' which are absent in the skeleton of 'egraph'
#' \item structural Hamming distance ('SHD') between 2 graphs is computed as TP+FP+the number of edges with an error in direction
#' \item TPR equals TP/(TP+FN)
#' \item FPR equals FP/(TN+FP) (TN stands for true negative edges)
#' \item FPRn equals FP/(TP+FN)
#' \item FDR equals FP/(TP+FP)
#' }
#'
#'@param egraph an object of class \code{\link[graph]{graphNEL}} (package `graph'), representing the graph which should be compared to a ground truth graph or an ajacency matrix corresponding to the graph
#'@param truegraph an object of class \code{\link[graph]{graphNEL}} (package `graph'), representing the ground truth graph or an ajacency matrix corresponding to this graph
#'@param cpdag logical, if TRUE (FALSE by default) both graphs are first converted to their respective equivalence class (CPDAG); this affects SHD calculation
#'@param rnd integer, rounding integer indicating the number of decimal places (round) when computing TPR, FPR, FPRn and FDR
#'@return a named numeric vector 8 elements: SHD, number of true positive edges (TP), number of false positive edges (FP), number of false negative edges (FN), true positive rate (TPR),
#'false positive rate (FPR), false positive rate normalized to the true number of edges (FPRn) and false discovery rate (FDR)
#'@examples
#'Asiascore<-scoreparameters("bde", Asia)
#'\dontrun{
#'eDAG<-orderMCMC(Asiascore)
#'compareDAGs(eDAG$DAG,Asiamat)
#'}
#'@export
compareDAGs<-function(egraph, truegraph, cpdag=FALSE, rnd=2) {
  
   if(is.matrix(egraph)) egraph<-m2graph(egraph)
   if(is.matrix(truegraph)) truegraph<-m2graph(truegraph)
   
   skeleton1<-graph2skeleton(egraph)
   n<-ncol(skeleton1)
   skeleton2<-graph2skeleton(truegraph)
   
  numedges1<-sum(skeleton1)
  numedges2<-sum(skeleton2)
  
  TN<-(n*n-n)/2-numedges2
  
  diff2<-skeleton2-skeleton1
  
  res<-vector()
  res[1]<-numedges1-sum(diff2<0) #TP
  res[2]<-sum(diff2<0) #FP
  res[3]<-sum(diff2>0)#FN
  res[4]<-round(res[1]/numedges2, rnd) #TPR
  res[5]<-round(res[2]/TN, rnd) #FPR
  res[6]<-round(res[2]/numedges2, rnd) #FPRn
  res[7]<-round(res[2]/(res[1]+res[2]),rnd) #FDR
  if(cpdag) {
    res[8]<-pcalg::shd(dag2cpdag(egraph), dag2cpdag(truegraph))
  } else res[8]<-pcalg::shd(egraph, truegraph)
  names(res)<-c("TP", "FP", "FN", "TPR", "FPR", "FPRn", "FDR", "SHD")
  return(res)
}

#'Comparing two DBNs
#'
#'This function compares one (estimated) DBN structure to another DBN (true DBN). Comparisons for initial and transitional structures are returned separately if \code{equalstruct} equals \code{TRUE}.
#'
#'@param eDBN an object of class \code{\link[graph]{graphNEL}} (or an ajacency matrix corresponding to this DBN), representing the DBN which should be compared to a ground truth DBN 
#'@param trueDBN an object of class \code{\link[graph]{graphNEL}} (or an ajacency matrix corresponding to this DBN), representing the ground truth DBN 
#'@param n.dynamic number of dynamic variables in one time slice of a DBN
#'@param n.static number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first n.static columns of the matrix
#'@param struct option used to determine if the initial or the transitional structure should be compared; accaptable values are init or trans
#'@return a vector of 5: SHD, number of true positive edges, number of false positive edges, number of false negative edges and true positive rate
#'@examples
#'testscore<-scoreparameters("bge", DBNdata, DBN=TRUE, 
#'dbnpar=list(samestruct=TRUE, slices=5, b=3))
#'\dontrun{
#'DBNfit<-iterativeMCMC(testscore, moveprobs=c(0.11,0.84,0.04,0.01))
#'compareDBNs(DBNfit$DAG,DBNmat, struct="trans", n.dynamic=12, n.static=3)
#'}
#'@export
compareDBNs<-function(eDBN, trueDBN, n.dynamic, n.static, struct=c("init","trans")) {
  
  if(length(struct)>1) {
    struct<-"trans"
    warning("parameter struct was not defined, 'trans' used by default")
  }

  n<-n.static+n.dynamic
  matsize<-n.static+2*n.dynamic
  
  if(!is.matrix(eDBN)) {
    adj1<-graph2m(eDBN)
  } else {
    adj1<-eDBN
  }
  
  if(!is.matrix(trueDBN)) {
    adj2<-graph2m(trueDBN)
  } else {
    adj2<-trueDBN
  }
  
  adj1<-adj1[1:matsize,1:matsize]
  adj2<-adj2[1:matsize,1:matsize]
  
  skeleton1<-graph2skeleton(eDBN)
  skeleton2<-graph2skeleton(trueDBN)
  
  if(struct=="trans") {
  
    skeleton1[,1:n]<-0
    skeleton2[,1:n]<-0
    
    adj1[,1:n]<-0
    adj2[,1:n]<-0
    
  } else if (struct=="init"){
    
    adj1<-adj1[1:n,1:n]
    adj2<-adj2[1:n,1:n]

    skeleton1<-skeleton1[1:n,1:n]
    skeleton2<-skeleton2[1:n,1:n]
    
  } else {
    stop("'struc' must be either 'init' or 'trans'!" )
  }
  eDAG<-m2graph(adj1)
  trueDAG<-m2graph(adj2)
  
  numedges1<-sum(skeleton1)
  numedges2<-sum(skeleton2)
  
  diff2<-skeleton2-skeleton1
  res<-vector()
  res["SHD"]<-pcalg::shd(eDAG, trueDAG)
  res["TP"]<-numedges1-sum(diff2<0) #TP
  res["FP"]<-sum(diff2<0) #FP
  res["FN"]<-sum(diff2>0)#FN
  res["TPR"]<-res["TP"]/numedges2 #TPR
  return(res)
}

#'Deriving interactions matrix
#'
#'This transforms a list of possible interactions between proteins downloaded from STRING database
#'into a matrix which can be used for blacklisting/penalization in BiDAG.
#'
#' @param curnames character vector with gene names which will be used in \code{BiDAG} learning function
#' @param int data frame, representing a interactions between genes/proteins downloaded from STRING (\url{https://string-db.org/}); two columns are necessary 'node1' and 'node2'
#' @param mapping (optional) data frame, representing a mapping between 'curnames' (gene names, usually the column names of 'data') and gene names used in interactions downloaded from STRING (\url{https://string-db.org/}); two columns are necessary 'queryItem' and 'preferredName'
#' @param type character, defines how interactions will be reflected in the output matrix; \code{int} will result in a matrix whose entries equal 1 if interaction is present in the list of interactions \code{int} and 0 otherwise; \code{blacklist} results in a matrix whose entries equal 0 when interaction is present in the list of interactions and 1 otherwise;
#' \code{pf} results in a matrix results in a matrix whose entries equal 1 is interaction is present in the list of interactions \code{int} and \code{pf} otherwise$ "int" by default
#' @param pf penalization factor for interactions, needed if \code{type}=pf 
#'@return square matrix whose entries  correspond to the list of interactions and parameter \code{type}
#'@examples
#'curnames<-colnames(kirp)
#'intmat<-string2mat(curnames, mapping, interactions, type="pf")
#'@export
string2mat<-function(curnames, int, mapping=NULL, type=c("int"), pf=2) {
  
  if(is.null(mapping)) {
    mapping<-cbind(curnames,curnames)
    colnames(mapping)<-c("queryItem","preferredName")
    mapping<-as.data.frame(mapping)
  }
  
  rownames(mapping)<-mapping$queryItem
  n<-length(curnames)
  aliases<-as.character(mapping[curnames,]$preferredName)
  nagenes<-which(is.na(aliases))
  if(length(nagenes)>0) {
    aliases[nagenes]<-curnames[nagenes]
    warning(paste(curnames[nagenes], "were not found in mapping; no interactions inferred"))
  }
  space<-matrix(0,nrow=n,ncol=n)
  rownames(space)<-aliases
  colnames(space)<-aliases
  for(i in 1:n) {
    curnode<-colnames(space)[i]
    nodz1<-intersect(int[which(int$node1==curnode),]$node2,aliases)
    nodz2<-intersect(int[which(int$node2==curnode),]$node1,aliases)
    if(length(nodz1)>0){
      space[curnode,nodz1]<-1
    }
    if(length(nodz2)>0){
      space[curnode,nodz2]<-1
    }
  }
  colnames(space)<-curnames
  rownames(space)<-curnames
  
  if(type=="int") space else if(type=="blacklist") 1*(!space) else (pf-1)*(!space)+1

}


#returns a matrix of a CPDAG corresponding to a given DAG
dagadj2cpadj<-function(adj) {
  g<-m2graph(adj)
  cpg<-pcalg::dag2cpdag(g)
  return(graph2m(cpg))
}

#returns a symmetric matrix of a skeleton corresponding to a given CPDAG
#UPPER TRIANGULAR VERSION!
adjacency2skeleton<-function(adj) {
  skel<-1*(adj|t(adj))
  skel<-ifelse(upper.tri(skel)==TRUE,skel,0)
  return(skel)
}

#returns a list of edges corresponding to an adjacency matrix
adjacency2edgel<-function(adj,nodes=NULL) {
  l<-ncol(adj)
  if (is.null(nodes)) {
    if (!all(is.character(colnames(adj)))) {
      V <- c(1:l)
      edL <- vector("list", length=l)
      names(edL) <- sapply(V,toString)
    } else {
      edL <- vector("list", length=l)
      names(edL) <- colnames(adj)
      V<-colnames(adj)
    }
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))}
  else {
    V <- nodes
    edL <- vector("list", length=l)
    names(edL) <- V
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))
  }
  return(edL)
  
}


#Deriving an adjacency matrix of the skeleton of a graph
#
#This function derives the skeleton matrix corresponding to a graph object
#
#g graph, object of class \code{\link[graph]{graphNEL}} (package `graph')
#returns a symmetric square matrix whose dimensions are the number of nodes in the graph \code{g},
# where element \code{[i,j]} equals \code{1} if there is a directed edge from node \code{i} to node \code{j},
#  or from node \code{j} to node \code{i}, in the graph \code{g}, and \code{0} otherwise
# examples 
#myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2)
#graph2skeleton.m(myDAG)
graph2skeleton<-function(g,upper=TRUE,outmat=TRUE) {
  if(!is.matrix(g)) {
    adj<-graph2m(g)
    colnames(adj)<-g@nodes
    rownames(adj)<-g@nodes
  } else {
    adj<-g
  }
  skel<-1*(adj|t(adj))
  if(upper) {
  skel<-ifelse(upper.tri(skel)==TRUE,skel,0)
  }
  if(outmat) {
    return(skel)
  } else {
    return(m2graph(skel))
  }
}


getRepr<-function(pdag,dag) {
  n<-ncol(pdag)
  for(i in 1:n) {
    for(j in 1:n) {
      if(pdag[i,j]==pdag[j,i] & pdag[i,j]==1) {
        pdag[i,j]<-dag[i,j]
        pdag[j,i]<-dag[j,i]
      }
    }
  }
  if(orderdag(pdag)!="error1") {
    return(pdag)
  } else {
    stop("not possible to resolve cycles!")
  }
}


orderdag<-function(adj) {
  n<-ncol(adj)
  allnodes<-c(1:n)
  curnodes<-c(1)
  order<-c()
  cntr<-1
  while(length(curnodes)<n & cntr<n) {
    npar<-apply(adj,2,sum)
    curnodes<-which(npar==0)
    order<-c(setdiff(curnodes,order),order)
    adj[curnodes,]<-0
    cntr<-cntr+1
  }
  
  if(sum(adj)==0) return(order)
  else return("error1")
  
}

