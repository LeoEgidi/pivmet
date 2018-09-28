#' Generate Data from a Gaussian Nested Mixture
#'
#' Simulate N observations from a nested Gaussian mixture model with k pre-specified components.
#'
#'
#' @param N Sample size, data dimension.
#' @param k Number of mixture components.
#' @param Mu Initial mean vector/matrix.
#' @param stdev Initial standard deviations (for univariate mixtures).
#' @param Sigma.p1 Covariance matrix for the first mixture level (for bivariate mixtures only).
#' @param Sigma.p2 Covariance matrix for the second mixture level (for bivariate mixture only).
#' @param W Mixture weights for the two levels,vector.
#' @return
#'
#' \item{ \code{y}}{Data values.}
#'
#' @details
#'
#' The functions allows to simulate values from a double (nested) univariate
#' Gaussian mixture:
#'
#' \deqn{
#' (Y_i|Z_i=j) sim \sum_{s=1}^{2} p_{js}\, \mathcal{N}(\mu_{j}, \sigma^{2}_{js}),
#' }
#'
#' or from a bivariate nested Gaussian mixture;
#'
#' \deqn{
#' (Y_i|Z_i=j) sim \sum_{s=1}^{2} p_{js}\, \mathcal{N}_{2}(\bm{\mu}_{j}, \Sigma_{s}).
#' }
#'
#' \code{Mu} is the mean input vector/matrix, \code{stdev} is a \eqn{k x 2}
#' matrix for the standard deviations of the \eqn{k} groups. \code{Sigma.p1} and
#' \code{Sigma.p2} are the covariances matrices for subgroups 1 and 2,
#' respectively. \code{W} is a vector of dimension 2 for the subgroups weights.
#' @examples
#'
#' # Bivariate mixture simulation with three components
#'
#' N  <- 2000
#' k  <- 3
#' M1 <- c(-45,8)
#' M2 <- c(45,.1)
#' M3 <- c(100,8)
#' Mu <- matrix(rbind(M1,M2,M3),c(k,2))
#' stdev    <- cbind(rep(1,k), rep(200,k))
#' Sigma.p1 <- matrix(c(stdev[1,1],0,0,stdev[1,1]),
#' nrow=2, ncol=2)
#' Sigma.p2 <- matrix(c(stdev[1,2],0,0,stdev[1,2]),
#'  nrow=2, ncol=2)
#' W   <- c(0.2,0.8)
#' sim <- piv_sim(N,k,Mu,stdev,Sigma.p1,Sigma.p2,W)
#' plot(sim$y, xlab="y[,1]", ylab="y[,2]")
#'
#' @export

piv_sim <- function(N,
                    k,
                    Mu,
                    stdev,
                    Sigma.p1 = matrix(c(1,0,0,1),2,2, byrow = TRUE),
                    Sigma.p2 = matrix(c(100,0,0,100),2,2, byrow = TRUE),
                    W = c(0.5, 0.5)){
  # Generation---------------

  if(missing(stdev)){
    stdev <- matrix(cbind(rep(1,k), rep(100,k)))
  }

  if (is.vector(Mu)){
    true.group <- sample(1:k,N,replace=TRUE,prob=rep(1/k,k))
    Spike <- array()
    matrixpi <- matrix(rep(W,k), nrow=k, ncol=2, byrow = T)
    sotto.gruppi <- matrix(0, nrow=k, ncol=N)
    y <- c()

    for (i in 1:N){
      sotto.gruppi[,i] <- sample(1:2, k, replace=T,
        prob=matrixpi[true.group[i],])
      y[i] <- rnorm(1, mean=Mu[true.group[i]],
                       sd=stdev[sotto.gruppi[true.group[i],i]])
    }

  }else{


    true.group <- sample(1:k,N,replace=TRUE,prob=rep(1/k,k))
    Spike <- array(c(Sigma.p1,Sigma.p2), dim=c(2,2,2))
    # Probability matrix of subgroups
    matrixpi <- matrix(rep(W,k), nrow=k, ncol=2, byrow = T)
    sotto.gruppi <- matrix(0, nrow=k, ncol=N)

  for (i in 1:N){
    sotto.gruppi[,i] <- sample(1:2,k,replace=T,
      prob=matrixpi[true.group[i],])
  }

  # Simulation of N units from a mixture of mixtures
  y <- matrix(NA,nrow=N,ncol=2)

  for (i in 1:length(true.group)){
    y[i,] <- mvrnorm(1, Mu[true.group[i],],
      Sigma=Spike[,,sotto.gruppi[true.group[i],i]])
  }
}

return(list(y=y, true.group=true.group, subgroups=sotto.gruppi))
}


