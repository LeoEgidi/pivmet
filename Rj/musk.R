#' K-means Clustering Using MUS algorithm
#'
#' Perform k-means clustering on a data matrix using MUS algorithm for seeding initialization.
#'
#' @param x A numeric matrix of data, or an object that can be coerced to such a matrix (such as a numeric vector or a dataframe with all numeric columns).
#' @param centers The number of clusters in the solution.
#' @param piv.criterion The pivotal criterion used for detecting pivotal units. If \code{centers <= 4}, default method is \code{MUS}. If \code{centers > 4}, the user may choose among the following: \code{maxsumint, maxsumnoint, maxsumdiff}.
#' @param iter.max The maximum number of iterations allowed.
#' @param num.seeds	The number of different starting random seeds to use. Each random seed results in a different k-means solution.
#' @param iter.mus The number of different ensembles for the MUS algorithm (if NULL, default is 1000)
#' @param prec.par The precision parameter used in the MUS algorithm
#' @param alg.type The type of clustering used for the initial seeding. Possible choices: \code{kmeans}, \code{KMeans}, \code{hclust}.
#'
#'@return A list with components
#'
#'\item{\code{cluster}}{A vector of integers indicating the cluster to which each point is allocated.}
#'\item{\code{centers}}{A matrix of cluster centres (centroids).}
#'\item{\code{totss}}{The total sum of squares.}
#'\item{\code{withinss}}{The within-cluster sum of squares for each cluster.}
#'\item{\code{tot.withinss}}{The within-cluster sum of squares summed across clusters.}
#'\item{\code{betwennss}}{The between-cluster sum of squared distances.}
#'\item{\code{size}}{ The number of points in each cluster.}
#'\item{\code{iter}}{The number of (outer) iterations.}
#'\item{\code{ifault}}{integer: indicator of a possible algorithm problem – for experts.}
#'\item{\code{pivots}}{The pivotal units identified by the MUS algorithm}
#'
#'@author Leonardo Egidi \url{legidi@units.it}
#'@examples
#'n  <- 620
#'k  <- 3
#'n1 <- 20
#'n2 <- 100
#'n3 <- 500
#'x  <- matrix(NA, n,2)
#'truegroup <- c( rep(1,n1), rep(2, n2), rep(3, n3))
#'
#'for (i in 1:n1){
#'  x[i,]=rmvnorm(1, c(1,5), sigma=diag(2))}
#'for (i in 1:n2){
#'  x[n1+i,]=rmvnorm(1, c(4,0), sigma=diag(2))}
#'for (i in 1:n3){
#'  x[n1+n2+i,]=rmvnorm(1, c(6,6), sigma=diag(2))}
#'
#'res <- MUSKMeans(x, k)
#'
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
#'    main="MUSK-means", cex.main=1.5)
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



MUSKMeans <- function(x, centers, piv.criterion,
  iter.mus, prec.par,
  alg.type, iter.max, num.seeds){
  #check on optional parameters
  if (missing(piv.criterion)){
    if (centers<=4 ){
      piv.criterion <- "MUS"
    }else{
      piv.criterion <- "maxsumdiff"
    }
  }
  if (missing(iter.mus)){
    iter.mus <- 1000
  }

  if (missing(alg.type)){
    alg.type <- "KMeans"
  }

  #type of clustering for initial clusters' assignment
  if (alg.type=="hclust"){
    cl <-cutree(hclust(dist(x), "average"),centers)
  }else if (alg.type=="KMeans"){
    cl <- KMeans(x,centers)$cluster
  }else if(alg.type=="kmeans"){
    cl <- kmeans(x,centers)$cluster
  }

  # tuning of precision MUS parameter
  if (missing(prec.par)){
    prec.par <- min( min(table(cl))-1, 5 )
  }

  #compute iter.mus different partitions
  if (is.vector(x)){
    n <- length(x)
  }else{
    n <- dim(x)[1]
  }
  H <- iter.mus
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
      prec.par <- prec.par
      mus_res  <- MUS(sim_matr, cl, prec.par)
      pivots   <- mus_res$maxima
    }else if (piv.criterion!="MUS"){

      #Other pivotal criteria
      z <- array(0,dim=c(n, centers, H))
      for (i in 1:H){
        for (j in 1:n){
          z[j,a[i,j],i] <- 1
        }
      }
      zm <- apply(z,c(1,3),FUN=function(x) sum(x*(1:length(x))))
      piv_sel <- pivotal_selection(Obj=c(1:7),
        k=centers, gIndex=cl,
        C=sim_matr, n=n, ZM=zm, maxima=c(1:centers),
        available_met = 7)
      if (piv.criterion=="maxsumint"){
        pivots <- piv_sel$Cg[,2]
      }else if(piv.criterion=="maxsumnoint"){
        pivots <- piv_sel$Cg[,5]
      }else if(piv.criterion=="maxsumdiff"){
        pivots <- piv_sel$Cg[,6]
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
    piv_sel <- pivotal_selection(Obj=c(1:6),
      k=centers, gIndex=cl,
      C=sim_matr, n=n, ZM=zm, maxima=pivots,
      available_met = 6)
    if (piv.criterion=="maxsumint"){
      pivots <- piv_sel$Cg[,2]
    }else if(piv.criterion=="maxsumnoint"){
      pivots <- piv_sel$Cg[,5]
    }else if(piv.criterion=="maxsumdiff"){
      pivots <- piv_sel$Cg[,6]
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
  d_mus   <- KMeans(x, centers=start)

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
