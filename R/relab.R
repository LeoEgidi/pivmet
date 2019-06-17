#' Performing the pivotal relabelling step and computing the relabelled posterior estimates
#'
#' This function allows to perform the pivotal relabelling procedure described in Egidi et al. (2018) and to obtain the relabelled posterior estimates.
#' @param mcmc The output of the MCMC sampling from \code{piv_MCMC}.
#'
#'@details
#'Prototypical models in which the label switching problem arises
#'are mixture models, as explained in the Details section of
#'the \code{piv_MCMC} function.
#'

#' These models are unidentified with respect to an arbitrary permutation
#' of the labels \eqn{1,...,k}. Relabelling means permuting
#'the labels at each iteration of the Markov chain in such
#'a way that the relabelled chain can be used to draw inferences
#'on component-specific parameters.
#'
#'
#' We assume here that a MCMC sample is obtained for the
#' posterior distribution of a Gaussian mixture model--for instance via
#' \code{piv_MCMC} function--with a prior distribution which is
#' labelling invariant.
#' Furthermore, suppose that we can find \eqn{k} units, one
#' for each group, which are (pairwise) separated with (posterior)
#' probability one
#' (that is, the posterior probability of any two of them being
#' in the same group
#' is zero).
#' It is then straightforward to use the \eqn{k} units,
#' called pivots in what follows and denoted by the indexes
#' \eqn{i_1,\ldots,i_k}, to identify the groups and to
#' relabel the chains:
#' for each MCMC iteration \eqn{h=1,\ldots, H} (\eqn{H} corresponds to
#' the argument \code{nMC}) and group
#'  \eqn{j=1,\ldots,k}, set
#'\deqn{
#'[\mu_j]_h=[\mu_{[Z_{i_{j}}]_h}]_h;
#'}
#'\deqn{
#'[Z_{i}]_h=j \mbox{ for } i:[Z_i]_h=[Z_{i_{j}}]_h.
#'}
#'The applicability of this strategy is limited by the existence of the pivots,
#'which is not guaranteed. The existence of the pivots is a requirement of the
#'method, meaning that its use is restricted to those chains—or
#'those parts of a chain—for which the pivots are present. First, although the
#'model is based on a mixture of \eqn{k} components, each iteration of the chain
#'may imply a different number of non-empty groups. Let then \eqn{[k]_h \leq k}
#'be the number of non-empty groups at iteration \eqn{h},
#'\deqn{
#'  [k]_h = \#\{j: [Z_i]_h=j\mbox{ for some }i\},
#'}
#'where \eqn{\#A} is the cardinality of the set \eqn{A}. Hence, the relabelling
#'procedure outlined above can be used only for the subset of the chain
#'for which \eqn{[k]_h=k}; let it be \deqn{\mathcal{H}_k=\{h:[k]_h= k\},}
#'which correspond to the argument \code{true.iter} given by \code{piv_MCMC}.
#'This means that the resulting relabelled chain is not a sample (of size \eqn{H})
#'from the posterior distribution, but a sample (of size \eqn{\#\mathcal{H}_k})
#'from the posterior
#'distribution conditional on there being (exactly) \eqn{k} non-empty groups.
#'Even if \eqn{k} non-empty groups are available, however,
#'there may not be \eqn{k} perfectly separated units. Let us define
#' \deqn{
#'  \mathcal{H}^{*}_k=\{ h\in\mathcal{H}_k : \exists r,s \mbox{ s.t. }
#'  [Z_{i_r}]_h=[Z_{i_s}]_h \}}
#'
#' that is, the set of iterations where (at least) two pivots are in the same
#' group.
#' In order for the pivot method to be applicable,
#' we need to exclude iterations \eqn{\mathcal{H}^{*}_k};
#' that is, we can perform the pivot relabelling on \eqn{\mathcal{H}_k-
#' \mathcal{H}^{*}_{k}}, corresponding to the argument \code{final_it}.
#'
#'

#' @return This function gives the relabelled posterior estimates--both mean and medians--obtained from the Markov chains of the MCMC sampling.
#'
#' \item{\code{final_it}}{The final number of valid MCMC iterations,
#' as explained in Details.}
#' \item{\code{final_it_p}}{The proportion of final valid MCMC iterations.}
#' \item{\code{rel_mean}}{The relabelled chains of the means: a \code{final_it}\eqn{\times k} matrix for univariate data,
#' or a \code{final_it}\eqn{\times 2 \times k} array for bivariate data.}
#' \item{\code{rel_sd}}{The relabelled chains of the sd's: a \code{final_it}\eqn{\times k} matrix for univariate data,
#' or a \code{final_it}\eqn{\times 2} matrix for bivariate data.}
#' \item{\code{rel_weight}}{The relabelled chains of the weights: a \code{final_it}\eqn{\times k} matrix.}
#'
#' @author Leonardo Egidi \url{legidi@units.it}
#' @references Egidi, L., Pappadà, R., Pauli, F. and Torelli, N. (2018). Relabelling in Bayesian Mixture
#'Models by Pivotal Units. Statistics and Computing, 28(4), 957-969.
#' @examples
#'
#' #Univariate simulation
#'\dontrun{
#' N   <- 250
#' nMC <- 2500
#' k   <- 3
#' p   <- rep(1/k,k)
#' x   <- 3
#' stdev <- cbind(rep(1,k), rep(20,k))
#' Mu    <- seq(-trunc(k/2)*x,trunc(k/2)*x,length=k)
#' W     <- c(0.2,0.8)
#' sim   <- piv_sim(N = N, k = k, Mu = Mu,
#'                  stdev = stdev, W=W)
#' res   <- piv_MCMC(y = sim$y, k =k, nMC = nMC)
#' rel   <- piv_rel(mcmc=res)
#'}
#'
#' #Bivariate simulation
#'\dontrun{
#' N <- 200
#' k <- 3
#' nMC <- 5000
#' M1  <- c(-.5,8)
#' M2  <- c(25.5,.1)
#' M3  <- c(49.5,8)
#' Mu  <- matrix(rbind(M1,M2,M3),c(k,2))
#' sds <- cbind(rep(1,k), rep(20,k))
#' Sigma.p1 <- matrix(c(sds[1,1]^2,0,0,sds[1,1]^2),
#'                    nrow=2, ncol=2)
#' Sigma.p2 <- matrix(c(sds[1,2]^2,0,0,sds[1,2]^2),
#'                    nrow=2, ncol=2)
#' W <- c(0.2,0.8)
#' sim <- piv_sim(N = N, k = k, Mu = Mu,
#'                Sigma.p1 = Sigma.p1,
#'                Sigma.p2 = Sigma.p2, W = W)
#' res <- piv_MCMC(y = sim$y, k = k, nMC = nMC)
#' rel <- piv_rel(mcmc = res)
#' piv_plot(y=sim$y, mcmc=res, rel_est = rel, type="chains")
#' piv_plot(y=sim$y, mcmc=res, rel_est = rel,
#'          type="hist")
#'}
#'
#' @export

piv_rel<-function(mcmc){

  ### checks


  ###


  N <- dim(mcmc$groupPost)[2]
  if (length(dim(mcmc$mcmc_mean_raw))==3){
  k <- dim(mcmc$mcmc_mean)[3]
  }else{
  k <- dim(mcmc$mcmc_mean_raw)[2]
  }
  nMC <- dim(mcmc$mcmc_mean_raw)[1]
  mu_switch <- mcmc$mcmc_mean
  tau_switch <- mcmc$mcmc_sd
  prob.st_switch <- mcmc$mcmc_weight
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
    tau_rel_median     <- c()  #vector of length k
    tau_rel_mean       <- c()
    tau_rel_median_tr  <- c()
    tau_rel_mean_tr    <- c()
    weights_rel_median     <- c()  #vector of length k
    weights_rel_mean       <- c()
    weights_rel_median_tr  <- c()
    weights_rel_mean_tr    <- c()
    groupD2           <- groupD[contD==0,]
    mu_switchD        <- mu_switch[contD==0,]
    tau_switchD       <- tau_switch[contD==0,]
    prob.st_switchD   <- prob.st_switch[contD==0,]
    true.iterD2       <- sum(contD==0)
    Final_It          <- true.iterD2/nMC
    mu_rel_complete   <- matrix(NA,true.iterD2, k)
    tau_rel_complete  <- matrix(NA, true.iterD2, k)
    weights_rel_complete <- matrix(NA, true.iterD2, k)


    if (true.iterD2!=0){
      for (m in 1:true.iterD2){
        for ( j in 1:k){
          vect_rel <- sort(mu_switchD[m,],
            decreasing=FALSE, index.return=TRUE)$ix
            mu_rel_complete[m,j] <-
              mu_switchD[m,
                vect_rel[groupD2[m, pivots[j]]] ]
            tau_rel_complete[m,j] <-
              tau_switchD[m,
                         vect_rel[groupD2[m, pivots[j]]] ]
            weights_rel_complete[m,j] <-
              prob.st_switchD[m,
                         vect_rel[groupD2[m, pivots[j]]] ]

          }
        }
      }else{
        stop("The number of MCMC iterations is too low, try increasing the argument nMC when you use the piv_MCMC function.")
        #mu_rel_median <- rep(NA,k)
        #mu_rel_mean   <- rep(NA,k)
      }

      mu_rel_median  <- apply(mu_rel_complete, 2, median)
      mu_rel_mean    <- apply(mu_rel_complete, 2, mean)
      mu_rel_median_tr  <- t(mu_rel_median)
      mu_rel_mean_tr    <- t(mu_rel_mean)
      tau_rel_median  <- apply(tau_rel_complete, 2, median)
      tau_rel_mean    <- apply(tau_rel_complete, 2, mean)
      tau_rel_median_tr  <- t(tau_rel_median)
      tau_rel_mean_tr    <- t(tau_rel_mean)
      weights_rel_median  <- apply(weights_rel_complete, 2, median)
      weights_rel_mean    <- apply(weights_rel_complete, 2, mean)
      weights_rel_median_tr  <- t(weights_rel_median)
      weights_rel_mean_tr    <- t(weights_rel_mean)

  }else{
    k <- dim(mu_switch)[3]
    mu_rel_median  <- array(NA,c(2,k))
    mu_rel_mean    <- array(NA,c(2,k))
    weights_rel_median <- c()
    weights_rel_mean <- c()
    groupD2        <- groupD[contD==0,]
    mu_switchD     <- mu_switch[contD==0,,]
    prob.st_switchD     <- prob.st_switch[contD==0,]
    true.iterD2    <- sum(contD==0)
    Final_It       <- true.iterD2/nMC
    mu_rel_complete  <- array(NA, dim=c(true.iterD2, 2,k))
    weights_rel_complete <- array(NA, dim=c(true.iterD2,k))

      if (true.iterD2!=0){
        for (m in 1:true.iterD2){
          for (j in 1:k){
            mu_rel_complete[m,,j] <-
              mu_switchD[m,,groupD2[m, pivots[j]]]
            weights_rel_complete[m,j] <-
              prob.st_switchD[m,groupD2[m, pivots[j]]]
          }
        }
      }else{
        stop("The number of MCMC iterations is too low, try increasing the argument nMC when you use the piv_MCMC function.")
        #mu_rel_median <- matrix(NA,c(2,k))
        #mu_rel_mean   <- matrix(NA,c(2,k))

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
  #weights_rel_complete[h,1] <- weights_rel_complete[h,ind[h,1,]]
  #weights_rel_complete[h,2] <- weights_rel_complete[h,ind[h,2,]]
  }
  mu_rel_median    <- apply(mu_rel_complete, c(2,3), median)
  mu_rel_mean      <- apply(mu_rel_complete, c(2,3), mean)
  mu_rel_median_tr <- t(mu_rel_median)
  mu_rel_mean_tr   <- t(mu_rel_mean)
  weights_rel_median    <- apply(weights_rel_complete, 2, median)
  weights_rel_mean      <- apply(weights_rel_complete, 2, mean)
  weights_rel_median_tr <- t(weights_rel_median)
  weights_rel_mean_tr   <- t(weights_rel_mean)
  tau_rel_median        <- apply(tau_switch, 2, median)
  tau_rel_complete      <- tau_switch
    }

 return(list( final_it = true.iterD2,
              rel_mean = mu_rel_complete,
              rel_sd = tau_rel_complete,
              rel_weight = weights_rel_complete,
              final_it_p = Final_It))
}

