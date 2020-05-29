#' Generate Data from a Gaussian Nested Mixture
#'
#' Simulate N observations from a nested Gaussian mixture model
#' with k pre-specified components under uniform group probabilities \eqn{1/k},
#' where each group is in turn
#' drawn from a further level consisting of two subgroups.
#'
#'
#' @param N The desired sample size.
#' @param k The desired number of mixture components.
#' @param Mu The input mean vector of length \eqn{k} for univariate
#' Gaussian mixtures; the input \eqn{k \times D} matrix with the
#' means' coordinates for multivariate Gaussian mixtures.
#' @param stdev For univariate mixtures, the  \eqn{k \times 2} matrix
#' of input standard deviations,
#' where the first column contains the parameters for subgroup 1,
#' and the second column contains the parameters for subgroup 2.
#' @param Sigma.p1 The \eqn{D \times D} covariance matrix for the first subgroup. For multivariate mixtures only.
#' @param Sigma.p2 The \eqn{D \times D} covariance matrix for the second subgroup. For multivariate mixtures only.
#' @param W The vector for the mixture weights of the two subgroups.
#' @return
#'
#' \item{\code{y}}{The \eqn{N} simulated observations.}
#' \item{\code{true.group}}{A vector of integers from \eqn{1:k}
#' indicating the values of the latent variables \eqn{Z_i}.}
#' \item{\code{subgroups}}{A \eqn{k \times N} matrix where
#' each row contains the index subgroup for the observations
#' in the \eqn{k}-th group.}
#'
#' @details
#'
#' The functions allows to simulate values from a double (nested) univariate
#' Gaussian mixture:
#'
#' \deqn{
#' (Y_i|Z_i=j) \sim \sum_{s=1}^{2} p_{js}\, \mathcal{N}(\mu_{j}, \sigma^{2}_{js}),
#' }
#'
#' or from a multivariate nested Gaussian mixture:
#'
#' \deqn{
#' (Y_i|Z_i=j) \sim \sum_{s=1}^{2} p_{js}\, \mathcal{N}_{D}(\bm{\mu}_{j}, \Sigma_{s}),
#' }
#'
#' where \eqn{\sigma^{2}_{js}} is the variance for the group \eqn{j} and
#' the subgroup \eqn{s} (\code{stdev} is the
#' argument for specifying the \code{k x 2} standard deviations
#' for univariate mixtures);
#'  \eqn{\Sigma_s} is the covariance matrix for the
#' subgroup \eqn{s, s=1,2}, where the two matrices are
#' specified by \code{Sigma.p1}
#' and \code{Sigma.p2} respectively; \eqn{\mu_j} and
#' \eqn{\bm{\mu}_j, \ j=1,\ldots,k}
#'   are the mean input vector and matrix respectively,
#'   specified by the argument \code{Mu};
#' \code{W} is a vector of dimension 2 for the subgroups weights.
#' @examples
#'
#' # Bivariate mixture simulation with three components
#'
#' N  <- 2000
#' k  <- 3
#' D <- 2
#' M1 <- c(-45,8)
#' M2 <- c(45,.1)
#' M3 <- c(100,8)
#' Mu <- matrix(rbind(M1,M2,M3),c(k,D))
#' sds <- cbind(rep(1,k), rep(20,k))
#' Sigma.p1 <- matrix(c( sds[1,1]^2, 0,0,
#'                       sds[1,1]^2), nrow=D, ncol=D)
#' Sigma.p2 <- matrix(c(sds[1,2]^2, 0,0,
#'                       sds[1,2]^2), nrow=D, ncol=D)
#' W   <- c(0.2,0.8)
#' sim <- piv_sim(N = N, k = k, Mu = Mu, Sigma.p1 = Sigma.p1,
#' Sigma.p2 = Sigma.p2, W = W)
#' graphics::plot(sim$y, xlab="y[,1]", ylab="y[,2]")
#'
#' @export

piv_sim <- function(N,
                    k,
                    Mu,
                    stdev,
                    Sigma.p1 = diag(2),
                    Sigma.p2 = 100*diag(2),
                    W = c(0.5, 0.5)){
  # Generation---------------


  #############
  ## checks

  # k

 if(is.numeric(Mu)==FALSE){
    stop("Specify a numeric vector or matrix for 'Mu'")
  }

  if (is.vector(Mu)){
    # check stdev
    if (k != length(Mu)){
      stop("The number of input means and the number of
      components do not match'")
    }

    if(missing(stdev)){
      stop("Argument 'stdev' missing with no default.")
      #stdev <- cbind(rep(1,k), rep(20,k))
    }

    if (dim(stdev)[1]!=k){
      stop("The number of rows of 'stdev' has to match
            the number of mixture components, k.")
    }
  }else if (is.matrix(Mu)){
    if (k != dim(Mu)[1]){
      stop("The number of input means and the number of
      components do not match")
    }

    # Sigma.p1 and Sigma.p2

    if (is.positive.definite(Sigma.p1)==FALSE |
        is.positive.definite(Sigma.p2)==FALSE){
      stop("Matrix covariances should be positive definite!")
    }

    # checks stdev
    if (missing(stdev)==FALSE){
      warning("'stdev' not required for bivariate data")
    }


  }

  # W

  if (sum(W)!=1){
    stop("Check that the sub-weights sum to one!")
  }else if(length(W)!=2){
    stop("The sub-weight vector should be of dimension two.")
  }


  ##########


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

    D <- dim(Mu)[2]
    if (dim(Sigma.p1)[2] + dim(Sigma.p2)[2] != 2*D ){
      stop("The number of dimensions in the covariance
           matrices is not well posed")
    }
    true.group <- sample(1:k,N,replace=TRUE,prob=rep(1/k,k))
    Spike <- array(c(Sigma.p1,Sigma.p2), dim=c(D,D,2))
    # Probability matrix of subgroups
    matrixpi <- matrix(rep(W,k), nrow=k, ncol=2, byrow = T)
    sotto.gruppi <- matrix(0, nrow=k, ncol=N)

  for (i in 1:N){
    sotto.gruppi[,i] <- sample(1:2,k,replace=T,
      prob=matrixpi[true.group[i],])
  }

  # Simulation of N units from a mixture of mixtures
  y <- matrix(NA,nrow=N,ncol=D)

  for (i in 1:length(true.group)){
    y[i,] <- mvrnorm(1, Mu[true.group[i],],
      Sigma=Spike[,,sotto.gruppi[true.group[i],i]])
  }
}

return(list(y=y, true.group=true.group, subgroups=sotto.gruppi))
}


