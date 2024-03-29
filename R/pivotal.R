#'Pivotal Selection via Co-Association Matrix
#'
#'
#'Finding pivotal units from a data partition and a
#'co-association matrix C
#'according to three different methods.
#'
#'@param C A \eqn{N \times N} co-association matrix, i.e.
#'a matrix whose elements are co-occurrences of pair of units
#'in the same cluster among \eqn{H} distinct partitions.
#'@param clusters A vector of integers from \eqn{1:k} indicating
#'a partition of the \eqn{N} units into, say, \eqn{k} groups.
#'
#'
#'
#'@details
#'
#' Given a set of \eqn{N} observations \eqn{(y_{1},y_{2},...,y_{N})}
#' (\eqn{y_i} may be a \eqn{d}-dimensional vector, \eqn{d \ge 1}),
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
#' within similarity (\code{"maxsumint"}) and the unit that
#' maximizes the difference between global within and
#' between similarities (\code{"maxsumdiff"}), respectively.
#' Alternatively, we may choose \eqn{i^{*} \in\mathcal J_j}, which minimizes:
#' \deqn{\sum_{p\not\in\mathcal J_j} c_{i^{*}p},}
#' obtaining the most distant unit among the members
#' that minimize the global dissimilarity between one group
#' and all the others (\code{"minsumnoint"}).
#' See the vignette for further details.
#'
#'@return
#'
#'\item{\code{pivots}}{A matrix with \eqn{k} rows and three
#' columns containing the indexes of the pivotal units for each method.}
#' @author Leonardo Egidi \email{legidi@units.it}
#' @references Egidi, L., Pappadà, R., Pauli, F. and Torelli, N. (2018). Relabelling in Bayesian Mixture
#'Models by Pivotal Units. Statistics and Computing, 28(4), 957-969.
#'
#'@examples
#' # Iris data
#'
#'data(iris)
#' # select the columns of variables
#'x<- iris[,1:4]
#'N <- nrow(x)
#'H <- 1000
#'a <- matrix(NA, H, N)
#'
#' # Perform H k-means partitions
#'
#'for (h in 1:H){
#'  a[h,] <- kmeans(x, centers = 3)$cluster
#'}
#' # Build the co-association matrix
#'
#'C <- matrix(NA, N,N)
#'for (i in 1:(N-1)){
#'  for (j in (i+1):N){
#'    C[i,j] <- sum(a[,i]==a[,j])/H
#'    C[j,i] <- C[i,j]
#'  }}
#'
#'km <- kmeans(x, centers =3)
#'
#' # Apply three pivotal criteria to the co-association matrix
#'
#' ris <- piv_sel(C, clusters = km$cluster)
#'
#' graphics::plot(iris[,1], iris[,2], xlab ="Sepal.Length", ylab= "Sepal.Width",
#' col = km$cluster)
#'
#'  # Add the pivots chosen by the maxsumdiff criterion
#'
#' points( x[ris$pivots[,3], 1:2], col = 1:3,
#' cex =2, pch = 8 )
#'
#'@export
#'

piv_sel<-function(C, clusters){

  ### checks

  # dimensions

  if (dim(C)[1]!=length(clusters)){
    stop("The length of the clusters vector does
         not coincide with the rows of the matrix C")
  }

  ###

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
        sum(C[j,com.gi],na.rm=TRUE),           #maxsumint
        sum(C[j,com.ngi]),                     #minsumnoint
        sum(C[j,com.gi],na.rm=TRUE)-sum(C[j,com.ngi])))
  }                                            #maxsumdiff

# Methods: from 1 to 4

    if (!is.null(ind.gi))
      Cg[g.i, 1] <- com.gi[which.max(ind.gi1[,2])]
    if (!is.null(ind.gi))
      Cg[g.i, 2] <- com.gi[which.min(ind.gi1[,3])]
    if (!is.null(ind.gi))
      Cg[g.i, 3] <- com.gi[which.max(ind.gi1[,4])]
    if (!is.null(ind.gi)){
      Cg1[g.i] <- ind.gi[do.call(order,
        as.data.frame(-ind.gi[,-1]))[1],1]
    }
  }

# For each method, we store the selected pivotal units
Cg <- Cg[,1:3]
colnames(Cg) <- c("maxsumint", "minsumnoint", "maxsumdiff")

  return(list(pivots=Cg))
}


