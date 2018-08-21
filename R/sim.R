#' Generate Data from a Nested Mixture
#'
#' Simulate N observations from a nested mixture model with k pre-specified components.
#'
#'
#' @param N sample size, data dimension.
#' @param k number of mixture components.
#' @param Mu initial mean vector/matrix.
#' @param stdev initial standard deviations (for univariate mixtures).
#' @param Sigma.p1 covariance matrix for the first mixture level.
#' @param Sigma.p2 covariance matrix for the second mixture level.
#' @param W mixture weights for the two levels,vector.
#' @return
#'
#' \code{y}  Data values.
#'
#' @examples
#'
#' Bivariate mixture simulation with three components
#'
#' N <- 200
#' k <- 3
#' M1 <- c(-.5,8)
#' M2 <- c(25.5,.1)
#' M3 <- c(49.5,8)
#' Mu <- matrix(rbind(M1,M2,M3),c(k,2))
#' stdev <- cbind(rep(1,k), rep(200,k))
#' Sigma.p1 <- matrix(c(stdev[1,1],0,0,stdev[1,1]),
#' nrow=2, ncol=2)
#' Sigma.p2 <- matrix(c(stdev[1,2],0,0,stdev[1,2]),
#'  nrow=2, ncol=2)
#' W <- c(0.2,0.8)
#' sim <- sim_mixture(N,k,Mu,stdev,Sigma.p1,Sigma.p2,W)
#'
#' @export

sim_mixture <- function(N,k,Mu, stdev,Sigma.p1,Sigma.p2,W){
  # Generation---------------

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


