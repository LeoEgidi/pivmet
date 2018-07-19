#' Perfrom the pivotal relabelling step and compute the relabelled posterior estimates
#'
#' This function allows to perform the pivotal relabelling procedure described in Egidi et al. (2017) and to obtain the relabelled posterior estimates.
#'
#' @param mu_switch The post-processed MCMC chains for the mean parameters.
#' @param group The units' group membership obtained after post-processing the chain.
#' @param clustering The output from clustering (\code{diana} or \code{HClust}).
#' @param Mu the true means (or an estimate for them).
#' @param nMC the number of total MCMC iterations (given in input for the function \code{BayesMCMC}).
#'
#'@details
#'Prototypical models in which the label switching problem arises are mixture models, where for a sample \eqn{\vy=(y_1,\ldots,y_n)} we assume
#'  \deqn{
#' (Y_i|Z_i=j) \sim f(y;\mu_j,\phi),
#' }
#' where the \eqn{Z_i}, \eqn{i=1,\ldots,n}, are i.i.d.\ random variables, \eqn{g=j,\dots,k},
#' \eqn{\phi} is a parameter which is common to all components,  \eqn{Z_i\in\{1,\ldots,k\}},
#' and
#'\deqn{
#' P(Z_i=k)=\pi_k.
#' }
#' we assume that an MCMC sample is obtained from the posterior distribution for model above--for instance via \code{bayesMCMC} function--with a prior distribution which is labelling invariant.
#'We denote as \eqn{\{\ca{\theta}:h=1,\ldots,H\}} the sample for the parameter \eqn{\theta}
#'=(\vmu,\vpi,\phi)\}, \eqn{H} being the number of MCMC iterations.
#'We assume that also a MCMC sample for the variable \eqn{Z} is obtained and denote it by  \eqn{\{\ca{Z}:h=1,\ldots,H\}}.
#' Let, then, \eqn{\mathcal J_1,\ldots,\mathcal J_{k}} be a partition, obtained via \code{diana} or \code{hclust}.
#' Furthermore, suppose that we can find \eqn{k} units, \eqn{i_1,\ldots,i_{k}}, one
#' for each group, which are (pairwise) separated with (posterior) probability one
#' (that is, the posterior probability of any two of them being in the same group
#' is zero). In terms of the matrix \eqn{C} with elements \eqn{c_{ip}=P(Z_i=Z_p|\mathcal y
#' )}, the \eqn{k\times k} submatrix with only the rows and columns correspondi
#' ng to \eqn{i_1,\ldots,i_{k}}, will be the identity matrix.
#' It is then straightforward to use the \eqn{k} units, called pivots in what follows, to identify the groups and to relabel the chains: for each \eqn{h=1,\ldots, H} and \eqn{j=1
#', \ldots,k}, set
#' \deqn{
#' \ca{\mu_j}=\ca{\mu_{\ca{Z_{i_j}}}};
#' }
#' \deqn{
#' \ca{Z_i}=j \mbox{ for } i:\ca{Z_i}=\ca{Z_{i_j}}.
#' }
#'First, although the model is based on a mixture of \eqn{k} components, each iteration of the chain may imply a different number of non-empty groups (that is, it may be that \eqn{[Z_i]_h\neq k\;\; \forall i} for some \eqn{j,h}); let then \eqn{\ca{k}\leq k} be the number of non-empty groups at iteration \eqn{h},
#'\deqn{
#'  \ca{k} = \#\{j: \ca{Z_i}=j\mbox{ for some }i\},
#'  }
#' where \eqn{\#A} is the cardinality of the set \eqn{A}.
#' If \eqn{\ca{k}<k} for some \eqn{h}, there cannot be \eqn{k} perfectly separated units, and so there cannot be \eqn{k} pivots.
#' Hence, the relabelling procedure outlined above can be used only for the subset of the chain for which \eqn{\ca{k}=k}; let it be \eqn{\hh_k=\{h:\ca{k}= k\}}, the \code{true.iter} argument.
#' Even if \eqn{k} non-empty groups are available, however, there may not be \eqn{k} perfectly separated units. Let us define
#' \deqn{
#'  \hh^{*}_k=\{ h\in\hh_k : \exists r,s \mbox{ s.t. } \ca{Z_{i_r}}=\ca{Z_{i_s}} \}
#'  }
#' that is, the set of iterations where (at least) two pivots are in the same group.
#'In order for the pivot method to be applicable, we need to exclude iterations #'\eqn{\hh^*_k}; that is, we can perform the pivot relabelling on \eqn{\hh_k-\hh_k^*}.
#' @return This function gives the relabelled posterior estimates--both mean and medians--obtained from the Markov chains of the MCMC sampling.
#'
#' \item{\code{mu_rel_mean}}{ Estimated posterior means}
#' \item{\code{mu_rel_median}}{ Estimated posterior medians}
#' \item{\code{mu_rel_complete}}{Complete relabelled chains}
#'
#' @author Leonardo Egidi \url{legidi@units.it}
#' @references Egidi, L., Pappada, R., Pauli, F. and Torelli, N. (2018). Relabelling in Bayesian Mixture
#'Models by Pivotal Units. Statistics and Computing, 28(4), 957-969, DOI 10.1007/s11222-017-  9774-2.
#' @examples
#'
#' #Univariate simulation
#'
#' N <- 250
#' nMC <- 2500
#' k <- 4
#' p <- rep(1/k,k)
#' x <- 3
#' stdev <- cbind(rep(1,k), rep(200,k))
#' Mu <- seq(-trunc(k/2)*x,trunc(k/2)*x,length=k)
#' W <- c(0.2,0.8)
#' sim <- sim_mixture(N,k,Mu,stdev,W=W)
#' output_bayes <- bayesMCMC(sim$y, k, nMC)
#' relab_est <-
#' pivotal_relabelling(mu_switch = output_bayes$mu_switch,
#'                     group = output_bayes$groupPost,
#'                     clustering= output_bayes$clust_sel,
#'                     Mu=output_bayes$Mu,
#'                     nMC = nMC)
#'
#'
#' plot_pivotal(y= sim$y,
#'              est = relab_est$mu_rel_median,
#'              chains=relab_est$mu_rel_complete,
#'              type="chains",
#'              mu_switch=output_bayes$mu_switch,
#'              n.iter=relab_est$Final_it,
#'              true.means= output_bayes$Mu)
#'
#' plot_pivotal(y= sim$y,
#'              est = relab_est$mu_rel_median,
#'              chains=relab_est$mu_rel_complete,
#'              type="estimates",
#'              mu_switch=output_bayes$mu_switch,
#'              n.iter=relab_est$Final_it,
#'              true.means= output_bayes$Mu)
#'
#' plot_pivotal(y= sim$y,
#'              est = relab_est$mu_rel_median,
#'              chains=relab_est$mu_rel_complete,
#'              type="estimates_hist",
#'              mu_switch=output_bayes$mu_switch,
#'              n.iter=relab_est$Final_it,
#'              true.means= output_bayes$Mu)
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
#' sim <- sim_mixture(N,k,Mu,stdev,Sigma.p1,Sigma.p2,W)
#' output_bayes <- bayesMCMC(sim$y, k, nMC)
#' relab_est <-
#' pivotal_relabelling(mu_switch=output_bayes$mu_switch,
#'                  group=output_bayes$groupPost,
#'                  clustering=output_bayes$clust_sel,
#'                  Mu=output_bayes$Mu,
#'                  nMC = nMC)
#'
#' plot_pivotal(y= sim$y,
#'              est = relab_est$mu_rel_median,
#'              chains=relab_est$mu_rel_complete,
#'              type="chains",
#'              mu_switch=output_bayes$mu_switch,
#'              n.iter=relab_est$Final_it,
#'              true.means= output_bayes$Mu)
#'
#' plot_pivotal(y= sim$y,
#'              est = relab_est$mu_rel_median,
#'              chains=relab_est$mu_rel_complete,
#'              type="estimates_hist",
#'              mu_switch=output_bayes$mu_switch,
#'              n.iter=relab_est$Final_it,
#'              true.means= output_bayes$Mu)
#'
#'
#'
#' @export

pivotal_relabelling<-function(mu_switch, group, clustering,
  Mu, nMC ){

  true.iter <- dim(mu_switch)[1]
  groupD <- array(NA, dim=c(true.iter, N))

  # cycle on number of iterations
  for (i in 1:true.iter){
      # cycle on number of groups
    for (j in 1:k){
     groupD[i, group[i,]==group[i, clustering$Cg[j]]] <- j
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
                vect_rel[groupD2[m, clustering$Cg[j]]] ]
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
              mu_switchD[m,,groupD2[m,clustering$Cg[j]]]
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

