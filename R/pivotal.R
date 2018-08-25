#'Pivotal Selection via Co-Association Matrix
#'
#'
#'Finding the pivots according to four different
#'methods involving a co-association matrix C. This is an internal function launched by \code{piv_MCMC}.
#'@param Obj Numerical string for the allowed pivotal criterion.
#'@param k The number of mixture components/groups.
#'@param gIndex Clusters' allocation.
#'@param C Co-association matrix.
#'@param n Data sample size
#'@param ZM Auxiliary matrix used for building \code{C}.
#'@param available_met Available criteria methods (integer).
#'@param maxima Initial assignment for MUS algorithm.
#'
#'@return
#'
#'\item{\code{Cg}}{ The pivotal units. }
#'
#'
#'@export
#'


piv_sel<-function(Obj, k, gIndex, C, n, ZM, maxima, available_met){

  if (missing(maxima)){
    maxima=c(1:k)
  }

Cg1 <- rep(NA, k)
Cg  <- matrix(NA, ncol=available_met, nrow=k)

 for (g.i in 1:k){
    com.gi  <-  (1:n)[gIndex==g.i]
    com.ngi <-  (1:n)[gIndex!=g.i]
    ind.gi  <- c()
    ind.gi1 <- c()


# choose the pivots within each group
  for (j in com.gi){
      ind.gi  <- rbind(ind.gi,c(j,sort(C[j,com.gi],
        decreasing=TRUE)))
      ind.gi1 <- rbind(ind.gi1,c(j,
        # criteria involving maximization
        #max(C[j,com.gi],na.rm=TRUE),          #dropped
        sum(C[j,com.gi],na.rm=TRUE),           #maxsumint
        # criteria involving minimization
        #min(C[j,com.gi],na.rm=TRUE),          #dropped
        #min(C[j,com.ngi]),                    #dropped
        sum(C[j,com.ngi]),                     #maxsumnoint
        # another criterion involving maximization
        sum(C[j,com.gi],na.rm=TRUE)-sum(C[j,com.ngi])))
  }                                            #maxsumdiff

# Methods: from 1 to 4
    # if (!is.null(ind.gi))
    #   Cg[g.i, 1] <- com.gi[which.max(ind.gi1[,2])]
    if (!is.null(ind.gi))
      Cg[g.i, 1] <- com.gi[which.max(ind.gi1[,2])]
    if (!is.null(ind.gi))
      Cg[g.i, 2] <- com.gi[which.min(ind.gi1[,3])]
    if (!is.null(ind.gi))
      Cg[g.i, 3] <- com.gi[which.min(ind.gi1[,4])]
    # if (!is.null(ind.gi))
    #   Cg[g.i, 5] <- com.gi[which.min(ind.gi1[,6])]
    # if (!is.null(ind.gi))
    #   Cg[g.i, 6] <- com.gi[which.max(ind.gi1[,7])]
    if (available_met==4){
    if (!is.null(ind.gi))
      Cg[, 4]    <-  t(maxima)
    }
    if (!is.null(ind.gi)){
      Cg1[g.i] <- ind.gi[do.call(order,
        as.data.frame(-ind.gi[,-1]))[1],1]
    }
  }

# For each method, we store the selected pivotal units
Cg <- Cg[,Obj]
# group1 contains the observation assigments to the groups obtained via pivots
group1 <- 0*ZM
# cycle on iterations
  for (i in 1:ncol(ZM)){
# cycle on number of groups
    for (j in 1:k){
      if (!is.na(Cg[j])){
        group1[ZM[,i] ==ZM[Cg[j],i],i] <- j
      }
    }
  }

# definition of the probabilities to belong to the groups for each unit
pr <- matrix(NA,nrow=k,ncol=n)
 for (kk in 1:k){
  pr[kk,] <- apply(group1,1,FUN=function(x) sum(x==kk)/length(x))
 }

# definition of the submatrix corresponding to the pivotal units

submatrix <- round(C[Cg,Cg],5)
T <- max(submatrix[upper.tri(submatrix)])

  return(list(pr=pr, Cg=Cg, Submatrix=submatrix, Max=T))
}


