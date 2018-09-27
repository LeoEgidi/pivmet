#'Pivotal Selection via Co-Association Matrix
#'
#'
#'Finding the pivots according to three different
#'methods involving a co-association matrix C.
#'@param C A \eqn{N x N} co-association matrix, i.e.
#'a matrix whose elements are co-occurences of pair of units
#'in the same cluster among \eqn{H} distinct partitions.
#'@param clusters A vector of integers indicating
#'a partition of the \eqn{N} units into, say, \eqn{k} groups.

#'
#'
#'@details
#'
#' Given a set of \eqn{N} observations \eqn{(y_{1},y_{2},...,y_{N})}
#' (\eqn{y_i} may be a \eqn{d}-dimensional vector, \eqn{d \in \mathbb{N}}),
#' consider clustering methods to obtain \eqn{H} distinct partitions
#' into \eqn{k} groups.
#' The matrix \code{C} is the co-association matrix,
#' where \eqn{c_{i,p}=n_{i,p}/H}, with \eqn{n_{i,p}} the number of times
#' the pair \eqn{(y_{i},y_{p})} is assigned to the same
#' cluster among the \eqn{H} partitions.
#'
#' Let \eqn{j} be the group containing units \eqn{\mathcal J_j},
#' the user may choose \eqn{{i^*}\in\mathcal J_j} that
#' maximizes one of the quantities:
#' \deqn{
#'  \sum_{p\in\mathcal J_j} c_{{i^*}p}}
#'
#'  or
#'  \deqn{\sum_{p\in\mathcal J_j} c_{{i^*}p} - \sum_{j\not\in\mathcal J_j} c_{{i^*}p}.
#' }
#'
#' These methods give the unit that maximizes the global
#' within similarity (\code{"maxsumint"} and the unit that
#' maximizes the difference between global within and
#' between similarities \code{"maxsumdiff"}, respectively.
#' Alternatively, we may choose \eqn{i^{*} \in\mathcal J_j}, which minimizes:
#' \deqn{\sum_{p\not\in\mathcal J_j} c_{i^{*}p},}
#' obtaining the most distant unit among the members
#' that minimize the global dissimilarity between one group
#' and all the others (\code{"maxsumnoint"}).
#' See the vignette for further details.
#'
#'@return
#'
#'\item{\code{pivots}}{ The pivotal units. }
#'
#'
#'@export
#'


piv_sel<-function(C, clusters){

N <- dim(C)[1]
k <- length(unique(clusters))

Cg1 <- rep(NA, k)
Cg  <- matrix(NA, ncol=3, nrow=k)

 for (g.i in 1:k){
    com.gi  <-  (1:N)[clusters==g.i]
    com.ngi <-  (1:N)[clusters!=g.i]
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
    # if (available_met==4){
    # if (!is.null(ind.gi))
    #   Cg[, 4]    <-  t(maxima)
    # }
    if (!is.null(ind.gi)){
      Cg1[g.i] <- ind.gi[do.call(order,
        as.data.frame(-ind.gi[,-1]))[1],1]
    }
  }

# For each method, we store the selected pivotal units
Cg <- Cg[,1:3]
# group1 contains the observation assigments to the groups obtained via pivots
# group1 <- 0*Z
# # cycle on iterations
#   for (i in 1:ncol(Z)){
# # cycle on number of groups
#     for (j in 1:k){
#       if (!is.na(Cg[j])){
#         group1[Z[,i] ==Z[Cg[j],i],i] <- j
#       }
#     }
#   }
#
# # definition of the probabilities to belong to the groups for each unit
# pr <- matrix(NA,nrow=k,ncol=N)
#  for (kk in 1:k){
#   pr[kk,] <- apply(group1,1,FUN=function(x) sum(x==kk)/length(x))
#  }
#
# # definition of the submatrix corresponding to the pivotal units
#
# submatrix <- round(C[Cg,Cg],5)
# T <- max(submatrix[upper.tri(submatrix)])

  return(list(
    #pr=pr,
     pivots=Cg
    #Submatrix=submatrix,
    #Max=T
    ))
}


