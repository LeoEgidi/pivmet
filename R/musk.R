#' k-means Clustering Using Pivotal Algorithms For Seeding
#'
#' Perform classical k-means clustering on a data matrix using pivots as
#' initial centers.
#'
#' @param x A \eqn{N \times D} data matrix, or an object that can be coerced to such a matrix (such as a numeric vector or a dataframe with all numeric columns).
#' @param centers The number of clusters in the solution.
#' @param piv.criterion The pivotal criterion used for identifying one pivot
#' for each group. Possible choices are: \code{"MUS", "maxsumint", "minsumnoint",
#' "maxsumdiff"}.
#' If \code{centers <= 4}, the default method is \code{"MUS"};
#' otherwise, the default method is \code{"maxsumdiff"} (see the details and
#' the vignette).
#' @param H If \code{"MUS"} is selected, this is the number of
#' distinct k-means partitions used for building a \eqn{N \times N}
#' co-association matrix.
#' @param alg.type The clustering algorithm for the initial partition of the
#' \eqn{N} units into the desired number of clusters.
#' Possible choices are \code{"KMeans"} (default) and \code{"hclust"}.
#' @param opt List of optional arguments to be passed to \code{MUS} or \code{KMeans}.
#'
#'
#' @details
#'
#' The function implements a modified version of k-means which aims at
#' improving the clustering solution starting from a careful seeding.
#' In particular, it performs a pivot-based initialization step
#' using pivotal methods to find the initial centers
#' for the clustering procedure. The starting point consists of multiple
#' runs of the classical k-means (which uses random seeds via \code{Kmeans}
#' function of the \code{RcmdrMisc} package)
#' with a fixed number of clusters
#' in order to build the co-association matrix of data units.
#'
#'
#' @return A list with components
#'
#'\item{\code{cluster}}{A vector of integers indicating the cluster to which each point is allocated.}
#'\item{\code{centers}}{A matrix of cluster centers (centroids).}
#'\item{\code{totss}}{The total sum of squares.}
#'\item{\code{withinss}}{The within-cluster sum of squares for each cluster.}
#'\item{\code{tot.withinss}}{The within-cluster sum of squares summed across clusters.}
#'\item{\code{betwennss}}{The between-cluster sum of squared distances.}
#'\item{\code{size}}{ The number of points in each cluster.}
#'\item{\code{iter}}{The number of (outer) iterations.}
#'\item{\code{ifault}}{integer: indicator of a possible algorithm problem – for experts.}
#'\item{\code{pivots}}{The pivotal units identified by the selected pivotal criterion.}
#'
#'
#'
#'
#'@author Leonardo Egidi \url{legidi@units.it}
#'@references
#'
#'Egidi, L., Pappadà, R., Pauli, F., Torelli, N. (2018).
#'K-means seeding via MUS algorithm. Conference Paper,
#'Book of Short Papers, SIS2018, ISBN: 9788891910233.
#'
#'@examples
#'
#' # Data generated from a mixture of three bivariate Gaussian distributions
#'\dontrun{
#'N  <- 620
#'k  <- 3
#'n1 <- 20
#'n2 <- 100
#'n3 <- 500
#'x  <- matrix(NA, N,2)
#'truegroup <- c( rep(1,n1), rep(2, n2), rep(3, n3))
#'
#'for (i in 1:n1){
#'  x[i,]=rmvnorm(1, c(1,5), sigma=diag(2))}
#'for (i in 1:n2){
#'  x[n1+i,]=rmvnorm(1, c(4,0), sigma=diag(2))}
#'for (i in 1:n3){
#'  x[n1+n2+i,]=rmvnorm(1, c(6,6), sigma=diag(2))}
#'
#' # Apply piv_KMeans with MUS as pivotal criterion
#'
#'res <- piv_KMeans(x, k)
#'
#'# Apply piv_KMeans with maxsumdiff as pivotal criterion
#'
#'res2 <- piv_KMeans(x, k, piv.criterion ="maxsumdiff")
#'
#' # Plot the data and the clustering solution
#'
#'par(mfrow=c(1,2), pty="s")
#'colors_cluster <- c("grey", "darkolivegreen3", "coral")
#'colors_centers <- c("black", "darkgreen", "firebrick")
#'plot(x, col = colors_cluster[truegroup],
#'    bg= colors_cluster[truegroup], pch=21, xlab="x[,1]",
#'    ylab="x[,2]", cex.lab=1.5,
#'    main="True data", cex.main=1.5)
#'
#'plot(x, col = colors_cluster[res$cluster],
#'    bg=colors_cluster[res$cluster], pch=21, xlab="x[,1]",
#'    ylab="x[,2]", cex.lab=1.5,
#'    main="piv_KMeans", cex.main=1.5)
#'points(x[res$pivots[1],1], x[res$pivots[1],2],
#'    pch=24, col=colors_centers[1],bg=colors_centers[1],
#'    cex=1.5)
#'points(x[res$pivots[2],1], x[res$pivots[2],2],
#'    pch=24,  col=colors_centers[2], bg=colors_centers[2],
#'    cex=1.5)
#'points(x[res$pivots[3],1], x[res$pivots[3],2],
#'    pch=24, col=colors_centers[3], bg=colors_centers[3],
#'    cex=1.5)
#'points(res$centers, col = colors_centers[1:k],
#'    pch = 8, cex = 2)
#'    }



piv_KMeans <- function(x,
                       centers,
                       alg.type = c("KMeans", "hclust"),
                       piv.criterion = c("MUS", "maxsumint", "minsumnoint", "maxsumdiff"),
                       H = 1000,
                       opt = list(iter.max =10, num.seeds =10, prec.par =5)){
  #check on optional parameters

  if(missing(opt)){
    iter.max  <- 10
    num.seeds <- 10
    prec.par  <- 5
  }else{
    iter.max  <-  c(10, list(opt)$iter.max)[which(c(is.null(list(opt)$iter.max), !is.null(list(opt)$iter.max)))]
    num.seeds <-  c(10, list(opt)$num.seeds)[which(c(is.null(list(opt)$num.seeds), !is.null(list(opt)$num.seeds)))]
    prec.par  <-  c(5, list(opt)$prec.par)[which(c(is.null(list(opt)$prec.par), !is.null(list(opt)$prec.par)))]
  }
  if (missing(piv.criterion)){
    if (centers<=4 ){
      piv.criterion <- "MUS"
    }else{
      piv.criterion <- "maxsumdiff"
    }
  }
  if (missing(H)){
    H <- 1000
  }

  if (missing(alg.type)){
    alg.type <- "KMeans"
  }

  #type of clustering for initial clusters' assignment
  if (alg.type=="hclust"){
    cl <-cutree(hclust(dist(x), "average"),centers)
  }else if (alg.type=="KMeans"){
    cl <- KMeans(x,centers)$cluster
  }

  # tuning of precision MUS parameter
  # if (missing(prec.par)){
  #   prec.par <- min( min(table(cl))-1, 5 )
  # }

  #compute H different partitions
  if (is.vector(x)){
    n <- length(x)
  }else{
    n <- dim(x)[1]
  }

  a <- matrix(NA, H, n)

  for (h in 1:H){
    a[h,] <- kmeans(x,centers)$cluster
  }
  sim_matr <- matrix(1, n,n)

  for (i in 1:(n-1)){
    for (j in (i+1):n){
      sim_matr[i,j] <- sum(a[,i]==a[,j])/H
      sim_matr[j,i] <- sim_matr[i,j]
    }
  }

  if (centers <=4){
    if (piv.criterion=="MUS"){

      #MUS algorithm
      #prec.par <- prec.par
      mus_res  <- MUS(sim_matr, cl)
      pivots   <- mus_res$pivots
    }else if (piv.criterion!="MUS"){

      #Other pivotal criteria
      z <- array(0,dim=c(n, centers, H))
      for (i in 1:H){
        for (j in 1:n){
          z[j,a[i,j],i] <- 1
        }
      }
      zm <- apply(z,c(1,3),FUN=function(x) sum(x*(1:length(x))))
      sel <- piv_sel(C=sim_matr, clusters=cl)
      if (piv.criterion=="maxsumint"){
        pivots <- sel$pivots[,1]
      }else if(piv.criterion=="minsumnoint"){
        pivots <- sel$pivots[,2]
      }else if(piv.criterion=="maxsumdiff"){
        pivots <- sel$pivots[,3]
      }
    }
  }else{
    z <- array(0,dim=c(n, centers, H))
    for (i in 1:H){
      for (j in 1:n){
        z[j,a[i,j],i] <- 1
      }
    }
    zm <- apply(z,c(1,3),FUN=function(x) sum(x*(1:length(x))))
    sel <- piv_sel(C=sim_matr,  clusters=cl)
    if (piv.criterion=="maxsumint"){
      pivots <- sel$pivots[,1]
    }else if(piv.criterion=="minsumnoint"){
      pivots <- sel$pivots[,2]
    }else if(piv.criterion=="maxsumdiff"){
      pivots <- sel$pivots[,3]
    }
  }

  #Initial seeding
  if(is.vector(x)){
    dim_x <- 1
    start <- c()
    for (k in 1:centers){
      start[k] <- as.double(x[pivots[k]])
    }
  }else if(is.matrix(x)){
    dim_x <- dim(x)[2]
    start <- matrix(NA, centers, dim_x )
    for (k in 1:centers){
      start[k,] <- as.double(x[pivots[k],])
    }
  }


  #MUSKmeans
  d_mus   <- KMeans(x, centers=start, iter.max = iter.max,
                    num.seeds = num.seeds)

  return(list(cluster=d_mus$cluster,
    centers=d_mus$centers,
    totss=d_mus$totss,
    withinss=d_mus$withinss,
    tot.withinss=d_mus$tot.withinss,
    betweenss=d_mus$betweenss,
    size=d_mus$size,
    iter=d_mus$iter,
    ifaults=d_mus$ifault,
    pivots = pivots,
    piv.criterion = piv.criterion))
}
