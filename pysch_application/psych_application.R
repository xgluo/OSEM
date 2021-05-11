########### Psychological Application ###########
# Required packages
library(BiDAG)
library(pcalg)
library(bnlearn)
library(corrplot)

# Data processing
# NOT RUN
# datRogers <- read.csv("OCDRogers.csv")   ## import data
# head(datRogers)
# dim(datRogers)
# datRogers[,"appetite"] <- datRogers[,"incappetite"] - datRogers[,"decappetite"]
# datRogers[,"appetite"] <- datRogers[,"appetite"] - min(datRogers[,"appetite"])
# datRogers[,"weight"] <- datRogers[,"weightgain"] - datRogers[,"weightloss"]
# datRogers[,"weight"] <- datRogers[,"weight"] - min(datRogers[,"weight"])
# datRogers[,"incappetite"] <- NULL
# datRogers[,"decappetite"] <- NULL
# datRogers[,"weightgain"] <- NULL
# datRogers[,"weightloss"] <- NULL
# symptoms <- colnames(datRogers)
# symptoms <- c(symptoms[1:12],"appetite","weight",symptoms[13:22])
# datRogers <- datRogers[,symptoms]
# save(datRogers, file = "OCDRogers.RData")

load("OCDRogers.RData")
n <- ncol(datRogers)
N <- nrow(datRogers)

###########
# McNally result
netdata <- as.data.frame(apply(datRogers, 2, as.numeric))  
set.seed(123)
bootnet <- boot.strength(netdata, R = 500, algorithm = "hc", algorithm.args = list(restart = 5, perturb = 10), debug = TRUE)  
avgnet <- averaged.network(bootnet, threshold = 0.85)
cpdag_mcnally <- dag2cpdag(amat(avgnet))
for (i in c(1:n)) {
  for (j in c(1:n)) {
    if (i != j) {
      temp <- bootnet$from == symptoms[i] & bootnet$to == symptoms[j]
      cpdag_mcnally[i,j] <- cpdag_mcnally[i,j] * bootnet[temp, "strength"] * bootnet[temp, "direction"]
    }
  }
}

###########
# Bootstrapping
# NOT RUN
# source("OCDRogers_boot.R")
# load("OCDRogers_boot.RData")
# res_cpdag <- res
# for (i in c(1:3)) {
#   for (j in c(1:500)) {
#     res_cpdag[,,i,j] <- dag2cpdag(res[,,i,j])
#   }
# }
# save(res_cpdag,file = "OCDRogers_boot_cpdag_all.RData")

load("OCDRogers_boot_cpdag_all.RData")
res_BDe <- res_cpdag[,,1,]
res_BGe <- res_cpdag[,,2,]
res_OSEM <- res_cpdag[,,3,]

###########
# point estimates

OSEMfit <- ordinalStructEM(datRogers,
                           usrpar = list(penType = "other",
                                         L = 5,
                                         lambda = 6))

BGE <- scoreparameters("bge", datRogers,bgepar = list(am = 0.05))
BGEfit <- iterativeMCMC(BGE)

BDE <- scoreparameters("bdecat",datRogers,bdecatpar = list(chi = 1.5))
BDEfit <- iterativeMCMC(BDE)

cpdag_OSEM <- dag2cpdag(OSEMfit$DAG)
cpdag_BGe <- dag2cpdag(BGEfit$DAG)
cpdag_BDe <- dag2cpdag(BDEfit$DAG)

dev.off()
plot(as(cpdag_BDe, "graphNEL"),main = "BDe")
dev.off()
plot(as(cpdag_BGe, "graphNEL"),main = "BGe")
dev.off()
plot(as(cpdag_OSEM, "graphNEL"),main = "OSEM")

###########

get_strength <- function (c, cboot) {
  n <- nrow(c)
  for (i in c(1:n)) {
    for (j in c(1:n)) {
      if (c[i,j] & i != j) {
        c[i,j] <- mean(apply(cboot,3,function (A) if ((A[i,j] + A[j,i]) == 1) {A[i,j]} else {A[i,j] / 2}))
      }
    }
  }
  return(c)
}
cpdag_OSEM <- get_strength(cpdag_OSEM, res_OSEM)
cpdag_BDe <- get_strength(cpdag_BDe, res_BDe)
cpdag_BGe <- get_strength(cpdag_BGe, res_BGe)


###########
# Generate heatmaps
par(mfcol = c(2,2))
corrplot(cpdag_OSEM, method = "shade", is.corr = FALSE,
         tl.col = "grey", cl.lim = c(0,1),
         mar = c(1,0,0,0)+0.5,addgrid.col = "lightgrey")
title(sub = "OSEM")
corrplot(cpdag_mcnally, method = "shade", is.corr = FALSE,
         tl.col = "grey", cl.lim = c(0,1),
         mar = c(1,0,0,0)+0.5,addgrid.col = "lightgrey")
title(sub = "McNally et al. (2017)")
corrplot(cpdag_BDe, method = "shade", is.corr = FALSE,
         tl.col = "grey", cl.lim = c(0,1),
         mar = c(1,0,0,0)+0.5,addgrid.col = "lightgrey")
title(sub = "BDe")
corrplot(cpdag_BGe, method = "shade", is.corr = FALSE,
         tl.col = "grey", cl.lim = c(0,1),
         mar = c(1,0,0,0)+0.5,addgrid.col = "lightgrey")
title(sub = "BGe")