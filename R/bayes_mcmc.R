#' JAGS Sampling for Gaussian Mixture Models and Clustering via Co-Association Matrix.
#'
#' Perform MCMC JAGS sampling for Gaussian mixture models, post-process the chains and apply a clustering technique to the MCMC sample. Pivotal units for each group are selected among four alternative criteria.
#' @param y N-dimensional data vector/matrix.
#' @param k Number of mixture components.
#' @param priors
#' @param nMC Number of MCMC iterations for the JAGS function execution.
#' @param piv.criterion The pivotal criterion used for identifying one pivot
#' for each group. Possible choices are: \code{"MUS", "maxsumint", "maxsumnoint",
#' "maxsumdiff"}.
#' If \code{k <= 4}, the default method is \code{"MUS"};
#' otherwise, the default method is \code{"maxsumdiff"} (see the details and
#' the vignette).
#' @param clustering The clustering technique adopted for partitioning the
#' \code{N} observations into \code{k} groups. Possible choices: \code{"diana"} (default),
#' \code{"hclust"}.
#'
#' @details
#' The function fits univariate and bivariate Bayesian Gaussian mixture models of the form
#' (here for univariate only):
#' \deqn{(Y_i|Z_i=j) \sim f(y;\mu_j,\phi_j),}
#' where the \eqn{Z_i}, \eqn{i=1,\ldots,N}, are i.i.d. random variables, \eqn{j=1,\dots,k},
#' \eqn{\phi_j} is the variance,  \eqn{Z_i \in {1,\ldots,k }}, and
#' \deqn{P(Z_i=j)=\pi_j.}
#' The likelihood of the model is then
#' \deqn{L(y;\mu,\pi,\phi) = \prod_{i=1}^n \sum_{j=1}^k \pi_j \mathcal{N}(\mu_j,\phi),}
#'with \eqn{\mu=(\mu_{1},\dots,\mu_{k})} component-specific parameters and \eqn{\pi=(\pi_{1},\dots,\pi_{k})} mixture weights. Let \eqn{\nu} denote a permutation of \eqn{{ 1,\ldots,k }}, and let \eqn{\nu(\mu)= (\mu_{\nu(1)},\ldots,} \eqn{ \mu_{\nu(k)})}, \eqn{ \nu(\pi)=(\pi_{\nu(1)},\ldots,\pi_{\nu(k)})} be the corresponding permutations of \eqn{\mu} and \eqn{\pi}. Denote by \eqn{V} the set of all the permutations of the indexes \eqn{{1,\ldots,k }}, the likelihood above is invariant under any permutation \eqn{\nu \in V}, that is
#' \deqn{
#' L(y;\mu,\pi,\phi) = L(y;\nu(\mu),\nu(\pi),\phi).}
#' As a consequence, the model is unidentified with respect to an arbitrary permutation of the labels.
#' When Bayesian inference for the model is performed, if the prior distribution \eqn{p_0(\mu,\pi,\phi)} is invariant under a permutation of the indices, then so is the posterior. That is, if \eqn{p_0(\mu,\pi,\phi) = p_0(\nu(\mu),\nu(\pi),\phi)}, then
#'\deqn{
#' p(\mu,\pi,\phi| y) \propto p_0(\mu,\pi,\phi)L(y;\mu,\pi,\phi)}
#' is multimodal with (at least) \eqn{k!} modes.
#'
#' Priors are chosen as weakly informative. For univariate mixtures,
#' the specification is the same as the function \code{BMMmodel} of the
#' \code{bayesmix} package:
#'
#'  \deqn{\mu_j \sim \mathcal{N}(0, 1/B0inv)}
#'  \deqn{\phi_j \sim \text{invGamma}(nu0Half, nu0S0Half)}
#'  \deqn{\pi \sim \text{Dirichlet}(1,\ldots,1)}
#'  \deqn{S0 \sim \text{Gamma}(g0Half, g0G0Half),}
#'
#'  with default values: \eqn{B0inv=0.1, nu0Half =10, S0=2,
#'  nu0S0Half= nu0Half*S0,
#'  g0Half = 5e-17, g0G0Half = 5e-33}, in accordance with the default
#'  specification
#'  \code{priors=list(kind = "independence", parameter = "priorsFish",
#'  hierarchical = "tau")} (see \code{bayesmix} for further details and choices).
#'
#'For bivariate mixtures, the prior specification is the following:
#'
#'\deqn{ \bm{\mu}_j  \sim \mathcal{N}_2(\bm{mu}_0, S2)}
#'\deqn{ 1/\Sigma \sim \text{Wishart(S3, 3)}}
#'\deqn{\pi \sim \text{Dirichlet}(1,\ldots,1),}
#'
#'where \eqn{S2} and \eqn{S3} are diagonal matrices
#'with diagonal elements (the variances)
#'equal to 1e+05. The user may specify other values for the hyperparameters
#'\eqn{\bm{\mu}_0, S2, S3} via \code{priors} argument in such a way:
#'
#'\code{priors =list(mu0 = c(1,1), S2 = matrix(c(0.002,0,0, 0.1),2,2, byrow=TRUE),
#'S3 = matrix(c(0.1,0,0,0.1), 2,2, byrow =TRUE))}.
#'
#'
#' The function performs JAGS sampling using the \code{bayesmix} package for univariate Gaussian mixtures, and the \code{runjags}
#' package for bivariate Gaussian mixtures. After MCMC sampling, this function
#' clusters the units in \code{k} groups,
#' calls the \code{piv_sel()} function and yields the
#' pivots obtained from one among four different
#' methods (the user may specify one among them via \code{piv.criterion}
#' argument):
#'  \code{"maxsumint"}, \code{"maxsumnoint"}, \code{"maxsumdiff"} and \code{"MUS"}
#'  (available only if \code{k < 5})
#' (see the vignette for thorough details).
#'
#' @return The function gives the MCMC output, the clustering solutions and the pivotal indexes. Here is a complete list of outputs.
#'
#' \item{\code{Freq}}{  \code{k x 2} matrix where: the first column
#' reports the number of units allocated to each group
#' as given by JAGS program; the second
#' column reports the same number of units as given by the
#' chains' post-processing.}
#' \item{\code{true.iter}}{ The number of MCMC iterations for which
#' the number of JAGS groups exactly coincides with the prespecified
#' number of groups \code{k}.}
#' \item{\code{z} }{  \code{N x k x true.iter} array with values: 1,
#' if the \eqn{i}-th unit belongs to the \eqn{j}-th group at
#' the \eqn{h}-th iteration; 0, otherwise.}
#' \item{\code{ris}}{  MCMC output matrix as provided by JAGS.}
#' \item{\code{groupPost}}{ \code{true.iter x N} matrix
#' with values from \code{1:k} indicating the post-processed group allocation
#' vector.}
#' \item{ \code{mu_switch}}{  If \code{y} is a vector, a \code{true.iter x k}
#' matrix with the post-processed MCMC chains for the mean parameters; if
#' \code{y} is a matrix, a \code{true.iter x 2 x k} array with
#' the post-processed MCMC chains for the mean parameters.}
#' \item{\code{mu_raw}}{ If \code{y} is a vector, a \code{nMC x k} matrix
#' with the raw MCMC chains for the mean parameters as given by JAGS; if
#' \code{y} is a matrix, a \code{nMC x 2 x k} array with the raw MCMC chains
#' for the mean parameters as given by JAGS.}
#' \item{\code{C}}{Co-association matrix constructed from the MCMC sample.}
#' \item{\code{grr}}{Group vector allocation as provided by
#' \code{"diana"} or \code{"hclust"}.}
#' \item{\code{pivots}}{ The pivotal units identified by the
#' selected pivotal criterion.}
#' \item{\code{piv.criterion}}{ Gives the pivotal criterion used for identifying
#' the pivots.}
#'
#'
#' @author Leonardo Egidi \url{legidi@units.it}
#' @references Egidi, L., Pappadà, R., Pauli, F. and Torelli, N. (2018). Relabelling in Bayesian Mixture
#'Models by Pivotal Units. Statistics and Computing, 28(4), 957-969.
#' @examples
#'
#' # Bivariate simulation
#'
#' N   <- 200
#' k   <- 4
#' nMC <- 1000
#' M1  <-c(-.5,8)
#' M2  <- c(25.5,.1)
#' M3  <- c(49.5,8)
#' M4  <- c(63.0,.1)
#' Mu  <- matrix(rbind(M1,M2,M3,M4),c(4,2))
#' stdev    <- cbind(rep(1,k), rep(200,k))
#' Sigma.p1 <- matrix(c(stdev[1,1],0,0,stdev[1,1]), nrow=2, ncol=2)
#' Sigma.p2 <- matrix(c(stdev[1,2],0,0,stdev[1,2]), nrow=2, ncol=2)
#' W <- c(0.2,0.8)
#' sim <- piv_sim(N,k,Mu, stdev, Sigma.p1,Sigma.p2,W)
#' res <- piv_MCMC(y = sim$y, k =k, nMC = nMC)
#' #changing priors
#' res2 <- piv_MCMC(y = sim$y,
#'                  priors = list (
#'                  mu0=c(1,1),
#'                  S2 = matrix(c(0.002,0,0, 0.1),2,2, byrow=TRUE),
#'                  S3 = matrix(c(0.1,0,0,0.1), 2,2, byrow =TRUE)),
#'                  k = k, nMC = nMC)
#'
#'
#'
#' # Fishery data (bayesmix package)
#'
#' data(fish)
#' y <- fish[,1]
#' k <- 5
#' nMC <- 5000
#' res <- piv_MCMC(y = y, k = k, nMC = nMC)
#' # changing priors
#' res2   <- piv_MCMC(y = y,
#'                    priors = list(kind = "condconjugate",
#'                    parameter = "priorsRaftery",
#'                    hierarchical = "tau"),  k =k, nMC = nMC)
#'
#' @export





piv_MCMC <- function(y,
                     k,
                     priors,
                     nMC,
                     piv.criterion = c("MUS", "maxsumint", "maxsumnoint", "maxsumdiff"),
                     clustering = c("diana", "hclust")){

  # Conditions about data dimension----------------

  if (is.vector(y)){
     if (missing(priors)){
       priors = list(kind = "independence", parameter = "priorsFish",
                     hierarchical = "tau")
     }
    N <- length(y)
    # JAGS code------------------------

    # Initial values
    mu_inits<- c()
    clust_inits <- kmeans(y, k)$cluster
    for (j in 1:k){
      mu_inits[j]<-mean(y[clust_inits==j])
    }
    # Data
    burn <- 1000

    # Model
    mod.mist.univ <- BMMmodel(y, k = k, initialValues = list(S0 = 2),
      priors = priors)
    control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
      burn.in = burn, n.iter = nMC, seed = 10)
    ogg.jags <- JAGSrun(y, model = mod.mist.univ, control = control)
    # Parameters' initialization

    J <- 3
    mcmc.pars <- array(data = NA, dim = c(nMC-length(1:burn), k, J))
    mcmc.pars[ , , 1] <- ogg.jags$results[-(1:burn), (N+k+1):(N+2*k)]
    mcmc.pars[ , , 2] <- ogg.jags$results[-(1:burn), (N+2*k+1):(N+3*k)]
    mcmc.pars[ , , 3] <- ogg.jags$results[-(1:burn), (N+1):(N+k)]

    mu_pre_switch_compl <-  mcmc.pars[ , , 1]
    tau_pre_switch_compl <-  mcmc.pars[ , , 2]
    prob.st_pre_switch_compl <-  mcmc.pars[ , , 3]

    mu <- mcmc.pars[,,1]
    tau <- mcmc.pars[,,2]
    prob.st <- mcmc.pars[,,3]
    group <-  ogg.jags$results[-(1:burn), 1:N] #gruppi
    FreqGruppiJags <- table(group)
    numeffettivogruppi <- apply(group,1,FUN = function(x) length(unique(x)))

    if (sum(numeffettivogruppi==k)==0){
      return(print("MCMC has not never been able to identify the required number of groups and the process has been interrupted"))
      #return(1)
    }

    ##saved in the output
    ris_prel <- ogg.jags$results[-(1:burn),]
    ris <- ris_prel[numeffettivogruppi==k,]
    true.iter <- nrow(ris)
    group <- ris[,1:N]

    group.orig <- group
    verigruppi <- as.double(names(table(group)))

    cont <- 0
    for (verogruppo in verigruppi){
      cont <- cont+1
      group.orig[group==verogruppo] <- cont          #aggiorna contatore pivot
    }
    cont                                           #qualche dubbio su sta parte

    k.orig <- k
    if (cont>1){
      k <- cont
    }
    mu <- mu[,verigruppi]
    tau <- tau[,verigruppi]
    prob.st <- prob.st[,verigruppi]

    M <- nrow(group)
    group <- group*0
    mu_switch <- array(rep(0, true.iter*k), dim=c(true.iter,k))
    z <- array(0,dim=c(N, k, true.iter))

    for (i in 1:true.iter){
      perm <- sample(1:k,k,replace=FALSE)
      for (j in 1:k){
        #post-processing
        group[i,group.orig[i,]==j] <- perm[j]
      }
      mu_switch[i,] <- mu[i,perm]
      tau[i,] <- tau[i,perm]
      prob.st[i,] <- prob.st[i,perm]
    }

    for (i in 1:true.iter){
      for (j in 1:N){
        z[j,group[i,j],i] <- 1
      }
    }


  }else if (is.matrix(y)){
    N <- dim(y)[1]

    # JAGS code------------------------

    # Initial values
    if (missing(priors)){
    mu0 <- as.vector(c(0,0))
    S2 <- matrix(c(1,0,0,1),nrow=2)/100000
    S3 <- matrix(c(1,0,0,1),nrow=2)/100000
    }else{
      mu0 <- priors$mu0
      S2 <- priors$S2
      S3 <- priors$S3
    }

    # Data
    dati.biv <- list(y = y, N = N, k = k, S2= S2, S3= S3, mu0=mu0,
      onesRepNclust = rep(1,k))

    # Model
    mod.mist.biv<-"model{
    # Likelihood:

    for (i in 1:N){
    yprev[i,1:2]<-y[i,1:2]
    y[i,1:2] ~ dmnorm(muOfClust[clust[i],],tauOfClust)
    clust[i] ~ dcat(pClust[1:k] )
    }

    # Prior:

    for (g in 1:k) {
    muOfClust[g,1:2] ~ dmnorm(mu0[],S2[,])}
    tauOfClust[1:2,1:2] ~ dwish(S3[,],3)
    Sigma[1:2,1:2] <- inverse(tauOfClust[,])
    pClust[1:k] ~ ddirch( onesRepNclust)
  }"


    # Parameters' initialization
    clust_inits <- KMeans(y, k)$cluster
    #cutree(hclust(dist(y), "average"),k)
    mu_inits <- matrix(0,k,2)
    for (j in 1:k){
      mu_inits[j,] <- cbind(mean(y[clust_inits==j,1]), mean(y[clust_inits==j,2]))
    }
    #Reorder mu_inits according to the x-coordinate
    mu_inits <-
      mu_inits[sort(mu_inits[,1], decreasing=FALSE, index.return=TRUE)$ix,]

    init1.biv <- dump.format(list(muOfClust=mu_inits,
      tauOfClust= matrix(c(15,0,0,15),ncol=2),
      pClust=rep(1/k,k), clust=clust_inits))
    moni.biv <- c("clust","muOfClust","tauOfClust","pClust")

    mod   <- mod.mist.biv
    dati  <- dati.biv
    init1 <- init1.biv
    moni  <- moni.biv

    # Jags execution
    ogg.jags <- run.jags(model=mod, data=dati, monitor=moni,
      inits=init1, n.chains=3,plots=FALSE, thin=1,
      sample=nMC, burnin=1000)
    # Extraction
    ris <- ogg.jags$mcmc[[1]]

    # Post- process of the chains----------------------
    group <- ris[,grep("clust[",colnames(ris),fixed=TRUE)]
    M <- nrow(group)
    H <- list()

    mu_pre_switch_compl <- array(rep(0, M*2*k), dim=c(M,2,k))
    for (i in 1:k){
      H[[i]] <- ris[,grep("muOfClust",colnames(ris),fixed=TRUE)][,c(i,i+k)]
    }
    for (i in 1:k){
      mu_pre_switch_compl[,,i] <- as.matrix(H[[i]])
    }
    # Discard iterations
    numeffettivogruppi <- apply(group,1,FUN = function(x) length(unique(x)))
    ris <- ris[numeffettivogruppi==k,]
    true.iter <- nrow(ris)

    if (sum(numeffettivogruppi==k)==0){
      return(print("MCMC has not never been able to identify the required number of groups and the process has been interrupted"))
      #return(1)
    }else{
      L<-list()
      mu_pre_switch <- array(rep(0, true.iter*2*k), dim=c(true.iter,2,k))
      for (i in 1:k){
        L[[i]] <- ris[,grep("muOfClust",colnames(ris),fixed=TRUE)][,c(i,i+k)]
      }
      for (i in 1:k){
        mu_pre_switch[,,i] <- as.matrix(L[[i]])
      }
    }

    group <- ris[,grep("clust[",colnames(ris),fixed=TRUE)]
    FreqGruppiJags <- table(group)
    tau <- ris[,grep("tauOfClust[",colnames(ris),fixed=TRUE)]
    prob.st <- ris[,grep("pClust[",colnames(ris),fixed=TRUE)]
    group.orig <- group
    verigruppi <- as.double(names(table(group)))
    prob.st <- prob.st[,verigruppi]

    mu_pre_switch <- mu_pre_switch[,,verigruppi]

    # Switching Post
    cont <- 0
    for (l in verigruppi){
      cont <- cont+1
      group.orig[group==l] <- cont
    }
    k.orig <- k
    if (cont > 1){
      k <- cont
    }
    mu_switch <- array(rep(0, true.iter*2*k), dim=c(true.iter,2,k))
    group <- group*0
    z <- array(0,dim=c(N, k, true.iter))

    for (i in 1:true.iter){
      perm <- sample(1:k,k,replace=FALSE)
      for (j in 1:k){
        #post-processing
        group[i,group.orig[i,]==j] <- perm[j]
      }
      mu_switch[i,,] <- mu_pre_switch[i,,perm]
      #tau[i,] <- tau[i,perm]
      prob.st[i,] <- prob.st[i,perm]
    }

    for (i in 1:true.iter){
      for (j in 1:N){
        z[j,group[i,j],i] <- 1
      }
    }

}

  FreqGruppiJagsPERM <- table(group)
  Freq <- cbind(FreqGruppiJags,FreqGruppiJagsPERM)
  colnames(Freq) <- c("JAGS raw groups", "JAGS post-processed groups")



  # Similarity matrix based on MCMC sampling------------------------
  nz <- dim(z)[1]
  M <- dim(z)[3]
  C <- matrix(1,nz,nz)
  zm <- apply(z,c(1,3),FUN=function(x) sum(x*(1:length(x))))

  for (i in 1:(nz-1)){
    for (j in (i+1):nz){
      C[i,j] <- sum(zm[i,]==zm[j,])/M
      C[j,i] <- C[i,j]
    }
  }
  matdissim <- 1-C
  diag(matdissim) <- 0

  # Clustering on dissimilarity matrix-------------

  if (missing(clustering)){
    #clustering <- "diana"
    gr  <- diana(matdissim,diss=TRUE)
    grr <- cutree(gr, k)
  }else if(clustering =="diana"){
    gr  <- diana(matdissim,diss=TRUE)
    grr <- cutree(gr, k)
  }else if(clustering == "hclust"){
    gr  <- hclust(as.dist(matdissim))
    grr <- cutree(gr, k)
  }

  available_met <- 3

  piv.criterion.choices <- c("maxsumint", "maxsumnoint",
    "maxsumdiff")

  if (missing(piv.criterion)){
    piv.criterion <- "maxsumdiff"
  }

  if (piv.criterion=="maxsumint"||
      piv.criterion=="maxsumnoint"||
      piv.criterion=="maxsumdiff" ){

    piv.index <- (1:3)[piv.criterion.choices==piv.criterion]
    piv.index.pivotal <- c(1,2,3)
    available_met <- 3
    x <- c(1:available_met)
    prec.par.1 <- min(min(table(grr))-1,5)
    clust  <-  piv_sel(C=C, clusters=as.vector(grr))

    pivots <- clust$pivots[,piv.index.pivotal[piv.index]]
  }else if(piv.criterion=="MUS"){
      if (k <=4 & sum(C==0)!=0){

          prec.par.1 <- min(min(table(grr))-1,5)
          mus_res    <- MUS(C, grr, prec.par.1)
          clust  <-  mus_res$pivots

  }else{

    print("maxsumdiff criterion instead of MUS has been adopted due to
          computational efficiency")
    clust  <-  piv_sel(C=C,  clusters=as.vector(grr))
    pivots <- clust$pivots[,3]
  }
}



  return(list( Freq=Freq, true.iter = true.iter, z=z, Mu = mu_inits,
    ris=ris, groupPost=group,
    mu_switch=mu_switch,
    mu_raw=mu_pre_switch_compl,
    C=C, grr=grr, pivots = pivots,
    piv.criterion = piv.criterion))
  }
