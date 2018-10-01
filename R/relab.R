#' Perfroming the pivotal relabelling step and computing the relabelled posterior estimates
#'
#' This function allows to perform the pivotal relabelling procedure described in Egidi et al. (2018) and to obtain the relabelled posterior estimates.
#' @param mcmc The output of the MCMC sampling from \code{piv_MCMC}.
#' @param nMC The number of total MCMC iterations (given in input to the \code{piv_MCMC} function, or any function suited for MCMC sampling).
#'
#'@details
#'Prototypical models in which the label switching problem arises are mixture models, where for a sample \eqn{y=(y_1,\ldots,y_n)} we assume
#'  \deqn{
#' (Y_i|Z_i=j) \sim f(y;\mu_j,\phi),
#' }
#' where the \eqn{Z_i}, \eqn{i=1,\ldots,n}, are i.i.d. random variables, \eqn{j=1,\dots,k},
#' \eqn{\phi} is a parameter which is common to all components,  \eqn{Z_i \in {1,\ldots,k}},
#' and
#'\deqn{
#' P(Z_i=j)=\pi_j.
#' }
#'
#' This model is unidentified with respect to an arbitrary permutation of the labels \eqn{1,...,k}. Relabelling means permuting
#'the labels at each iteration of the Markov chain in such
#'a way that the relabelled chain can be used to draw inferences
#'on component-specific parameters.
#'
#'
#' We assume here that an MCMC sample is obtained from the posterior distribution for model above--for instance via \code{piv_MCMC} function--with a prior distribution which is labelling invariant.
#' Furthermore, suppose that we can find \eqn{k} units, one
#' for each group, which are (pairwise) separated with (posterior) probability one
#' (that is, the posterior probability of any two of them being in the same group
#' is zero).
#' It is then straightforward to use the \eqn{k} units, called pivots in what follows, to identify the groups and to relabel the chains (see the vignette for thorough details).

#' @return This function gives the relabelled posterior estimates--both mean and medians--obtained from the Markov chains of the MCMC sampling.
#'
#' \item{\code{mu_rel_mean}}{ \code{k}-vector or \code{k x 2}
#' matrix of estimated posterior means for the mean parameters.}
#' \item{\code{mu_rel_median}}{ \code{k}-vector or \code{k x 2}
#' matrix of estimated posterior medins for the mean parameters.}
#' \item{\code{mu_rel_complete}}{Complete relabelled chains}
#' \item{\code{Final_It}}{The final number of valid iterations}
#'
#' @author Leonardo Egidi \url{legidi@units.it}
#' @references Egidi, L., Pappada, R., Pauli, F. and Torelli, N. (2018). Relabelling in Bayesian Mixture
#'Models by Pivotal Units. Statistics and Computing, 28(4), 957-969, DOI 10.1007/s11222-017-  9774-2.
#' @examples
#'
#' #Univariate simulation
#'
#' N   <- 250
#' nMC <- 2500
#' k   <- 3
#' p   <- rep(1/k,k)
#' x   <- 3
#' stdev <- cbind(rep(1,k), rep(200,k))
#' Mu    <- seq(-trunc(k/2)*x,trunc(k/2)*x,length=k)
#' W     <- c(0.2,0.8)
#' sim   <- piv_sim(N,k,Mu,stdev,W=W)
#' res   <- piv_MCMC(sim$y, k, nMC)
#' rel   <- piv_rel(mcmc=res, nMC = nMC)
#'
#'
#' #Bivariate simulation
#'
#' N <- 200
#' k <- 3
#' nMC <- 5000
#' M1  <- c(-.5,8)
#' M2  <- c(25.5,.1)
#' M3  <- c(49.5,8)
#' Mu  <- matrix(rbind(M1,M2,M3),c(k,2))
#' stdev <- cbind(rep(1,k), rep(200,k))
#' Sigma.p1 <- matrix(c(stdev[1,1],0,0,stdev[1,1]),
#'                    nrow=2, ncol=2)
#' Sigma.p2 <- matrix(c(stdev[1,2],0,0,stdev[1,2]),
#'                    nrow=2, ncol=2)
#' W <- c(0.2,0.8)
#' sim <- piv_sim(N,k,Mu,stdev,Sigma.p1,Sigma.p2,W)
#' res <- piv_MCMC(sim$y, k, nMC)
#' rel <- piv_rel(mcmc = res, nMC = nMC)
#' piv_plot(y=sim$y, mcmc=res, rel_est = rel, type="chains")
#' piv_plot(y=sim$y, mcmc=res, rel_est = rel,
#'          type="estimates_hist")
#'
#'
#' @export

piv_rel<-function(mcmc, nMC ){

  mu_switch <- mcmc$mu_switch
  group <-  mcmc$groupPost
  pivots <- mcmc$pivots
  Mu <- mcmc$Mu

  true.iter <- dim(mu_switch)[1]
  groupD <- array(NA, dim=c(true.iter, N))

  # cycle on number of iterations
  for (i in 1:true.iter){
      # cycle on number of groups
    for (j in 1:k){
     groupD[i, group[i,]==group[i, pivots[j]]] <- j
        #group[i, clustering[[u]]$Cg[j]]
     }
    }
  # Finding final number of iterations : H_{G}-H^{*}_G, as explained     in the paper
  contD <- c()
  for (i in 1:true.iter){
      contD[i] <- sum(is.na(groupD[i,])==TRUE)
    }

# Final_It contains the final valid number of iterations


 if (length(dim(mu_switch))==2){
    k <- dim(mu_switch)[2]
    mu_rel_median     <- c()  #vector of length k
    mu_rel_mean       <- c()
    mu_rel_median_tr  <- c()
    mu_rel_mean_tr    <- c()
    groupD2           <- groupD[contD==0,]
    mu_switchD        <- mu_switch[contD==0,]
    true.iterD2       <- sum(contD==0)
    Final_It          <- true.iterD2/nMC
    mu_rel_complete   <- matrix(NA,true.iterD2, k)


    if (true.iterD2!=0){
      for (m in 1:true.iterD2){
        for ( j in 1:k){
          vect_rel <- sort(mu_switchD[m,],
            decreasing=FALSE, index.return=TRUE)$ix
            mu_rel_complete[m,j] <-
              mu_switchD[m,
                vect_rel[groupD2[m, pivots[j]]] ]
          }
        }
      }else{
        mu_rel_median <- rep(NA,k)
        mu_rel_mean   <- rep(NA,k)
      }

      mu_rel_median  <- apply(mu_rel_complete, 2, median)
      mu_rel_mean    <- apply(mu_rel_complete, 2, mean)
      mu_rel_median_tr  <- t(mu_rel_median)
      mu_rel_mean_tr    <- t(mu_rel_mean)

  }else{
    k <- dim(mu_switch)[3]
    mu_rel_median  <- array(NA,c(2,k))
    mu_rel_mean    <- array(NA,c(2,k))
    groupD2        <- groupD[contD==0,]
    mu_switchD     <- mu_switch[contD==0,,]
    true.iterD2    <- sum(contD==0)
    Final_It       <- true.iterD2/nMC
    mu_rel_complete  <- array(NA, dim=c(true.iterD2, 2,k))


      if (true.iterD2!=0){
        for (m in 1:true.iterD2){
          for (j in 1:k){
            mu_rel_complete[m,,j] <-
              mu_switchD[m,,groupD2[m, pivots[j]]]
          }
        }
      }else{
        mu_rel_median <- matrix(NA,c(2,k))
        mu_rel_mean   <- matrix(NA,c(2,k))

      }

    ind <- array(NA, c( nrow(mu_rel_complete) ,2, k))
      for (g in 1:nrow(mu_rel_complete)){
        prel1 <- c()
        prel2 <- c()
        for (h in 1:k){
    prel1[h] <- which.min((mu_rel_complete[g,1,]-Mu[h,1])^2)
    prel2[h] <- which.min((mu_rel_complete[g,2,]-Mu[h,2])^2)
         }
        ind[g,1,] <- prel1
        ind[g,2,] <- prel2
      }

    mu_rel_median_tr  <- array(NA, c(k,2))
    mu_rel_mean_tr    <- array(NA, c(k,2))

 for (h in 1:nrow(mu_rel_complete)){
  mu_rel_complete[h,1,] <- mu_rel_complete[h,1, ind[h,1,]]
  mu_rel_complete[h,2,] <- mu_rel_complete[h,2, ind[h,2,]]
      }
  mu_rel_median    <- apply(mu_rel_complete, c(2,3), median)
  mu_rel_mean      <- apply(mu_rel_complete, c(2,3), mean)
  mu_rel_median_tr <- t(mu_rel_median)
  mu_rel_mean_tr   <- t(mu_rel_mean)
    }

 return(list(mu_rel_median = mu_rel_median_tr,
              mu_rel_mean = mu_rel_mean_tr,
              mu_rel_complete = mu_rel_complete,
              Final_It = Final_It))
}

