#' JAGS/Stan Sampling for Gaussian Mixture Models and Clustering via Co-Association Matrix.
#'
#' Perform MCMC JAGS sampling or HMC Stan sampling for Gaussian mixture models, post-process the chains and apply a clustering technique to the MCMC sample. Pivotal units for each group are selected among four alternative criteria.
#' @param y \eqn{N}-dimensional vector for univariate data or
#' \eqn{N \times D} matrix for multivariate data.
#' @param k Number of mixture components.
#' @param nMC Number of MCMC iterations for the JAGS/Stan function execution.
#' @param priors Input prior hyperparameters (see Details for default options).
#' @param piv.criterion The pivotal criterion used for identifying one pivot
#' for each group. Possible choices are: \code{"MUS", "maxsumint", "minsumnoint",
#' "maxsumdiff"}.
#' The default method is \code{"maxsumint"} (see the Details and
#' the vignette).
#' @param clustering The algorithm adopted for partitioning the
#' \eqn{N} observations into \code{k} groups. Possible choices are \code{"diana"} (default) or
#' \code{"hclust"} for divisive and agglomerative hierarchical clustering, respectively.
#' @param software The selected MCMC method to fit the model: \code{"rjags"} for the JAGS method, \code{"rstan"} for the Stan method.
#' Default is \code{"rjags"}.
#' @param burn The burn-in period (only if method \code{"rjags"} is selected). Default is \code{0.5}\eqn{\times} \code{nMC}.
#' @param chains A positive integer specifying the number of Markov chains. The default is 4.
#' @param cores The number of cores to use when executing the Markov chains in parallel (only if
#' \code{software="rstan"}). Default is 1.
#' @param sparsity Allows for sparse finite mixtures, default is \code{TRUE}.
#'
#' @details
#' The function fits univariate and multivariate Bayesian Gaussian mixture models of the form
#' (here for univariate only):
#' \deqn{(Y_i|Z_i=j) \sim \mathcal{N}(\mu_j,\sigma_j),}
#' where the \eqn{Z_i}, \eqn{i=1,\ldots,N}, are i.i.d. random variables, \eqn{j=1,\dots,k},
#' \eqn{\sigma_j} is the group variance,  \eqn{Z_i \in {1,\ldots,k }} are the
#' latent group allocation, and
#' \deqn{P(Z_i=j)=\eta_j.}
#' The likelihood of the model is then
#' \deqn{L(y;\mu,\eta,\sigma) = \prod_{i=1}^N \sum_{j=1}^k \eta_j \mathcal{N}(\mu_j,\sigma_j),}
#' where \eqn{(\mu, \sigma)=(\mu_{1},\dots,\mu_{k},\sigma_{1},\ldots,\sigma_{k})}
#' are the component-specific parameters and \eqn{\eta=(\eta_{1},\dots,\eta_{k})}
#' the mixture weights. Let \eqn{\nu} denote a permutation of \eqn{{ 1,\ldots,k }},
#' and let \eqn{\nu(\mu)= (\mu_{\nu(1)},\ldots,} \eqn{ \mu_{\nu(k)})},
#' \eqn{\nu(\sigma)= (\sigma_{\nu(1)},\ldots,} \eqn{ \sigma_{\nu(k)})},
#' \eqn{ \nu(\eta)=(\eta_{\nu(1)},\ldots,\eta_{\nu(k)})} be the
#' corresponding permutations of \eqn{\mu}, \eqn{\sigma} and \eqn{\eta}.
#'  Denote by \eqn{V} the set of all the permutations of the indexes
#'  \eqn{{1,\ldots,k }}, the likelihood above is invariant under any
#'  permutation \eqn{\nu \in V}, that is
#' \deqn{
#' L(y;\mu,\eta,\sigma) = L(y;\nu(\mu),\nu(\eta),\nu(\sigma)).}
#' As a consequence, the model is unidentified with respect to an
#' arbitrary permutation of the labels.
#' When Bayesian inference for the model is performed,
#' if the prior distribution \eqn{p_0(\mu,\eta,\sigma)} is invariant under a permutation of the indices, then so is the posterior. That is, if \eqn{p_0(\mu,\eta,\sigma) = p_0(\nu(\mu),\nu(\eta),\sigma)}, then
#'\deqn{
#' p(\mu,\eta,\sigma| y) \propto p_0(\mu,\eta,\sigma)L(y;\mu,\eta,\sigma)}
#' is multimodal with (at least) \eqn{k!} modes.
#'
#' Depending on the selected software, the model parametrization
#' changes in terms of the prior choices.
#' Precisely, the JAGS philosophy with the underlying Gibbs sampling
#' is to use noninformative priors, and conjugate priors are
#' preferred for computational speed.
#' Conversely, Stan adopts weakly informative priors,
#' with no need to explicitly use the conjugacy.
#' For univariate mixtures, when
#' \code{software="rjags"} the specification is the same as the function \code{BMMmodel} of the
#' \code{bayesmix} package:
#'
#'  \deqn{\mu_j \sim \mathcal{N}(\mu_0, 1/B0inv)}
#'  \deqn{\sigma_j \sim \mbox{invGamma}(nu0Half, nu0S0Half)}
#'  \deqn{\eta \sim \mbox{Dirichlet}(1,\ldots,1)}
#'  \deqn{S0 \sim \mbox{Gamma}(g0Half, g0G0Half),}
#'
#'  with default values: \eqn{\mu_0=0, B0inv=0.1, nu0Half =10, S0=2,
#'  nu0S0Half= nu0Half\times S0,
#'  g0Half = 5e-17, g0G0Half = 5e-33}, in accordance with the default
#'  specification:
#'
#'  \code{priors=list(kind = "independence", parameter = "priorsFish",
#'  hierarchical = "tau")}
#'
#'  (see \code{bayesmix} for further details and choices).
#'
#'  When \code{software="rstan"}, the prior specification is:
#'
#'  \deqn{\mu_j \sim \mathcal{N}(\mu_0, 1/B0inv)}
#'  \deqn{\sigma_j \sim \mbox{Lognormal}(\mu_{\sigma}, \tau_{\sigma})}
#'  \deqn{\eta_j \sim \mbox{Uniform}(0,1),}
#'
#'  with default values: \eqn{\mu_0=0, B0inv=0.1, \mu_{\sigma}=0, \tau_{\sigma}=2}.
#' The users may specify new hyperparameter values with the argument:
#'
#' \code{priors=list(mu_0=1, B0inv=0.2, mu_sigma=3, tau_sigma=5)}
#'
#'For multivariate mixtures, when \code{software="rjags"} the prior specification is the following:
#'
#'\deqn{ \bm{\mu}_j  \sim \mathcal{N}_D(\bm{\mu}_0, S2)}
#'\deqn{ \Sigma^{-1} \sim \mbox{Wishart}(S3, D+1)}
#'\deqn{\eta \sim \mbox{Dirichlet}(\bm{\alpha}),}
#'
#'where  \eqn{\bm{\alpha}} is a \eqn{k}-dimensional vector
#'and \eqn{S_2} and \eqn{S_3}
#'are positive definite matrices. By default, \eqn{\bm{\mu}_0=\bm{0}},
#'\eqn{\bm{\alpha}=(1,\ldots,1)} and \eqn{S_2} and \eqn{S_3} are diagonal matrices,
#'with diagonal elements
#'equal to 1e+05. The user may specify other values for the hyperparameters
#'\eqn{\bm{\mu}_0, S_2, S_3} and \eqn{\bm{\alpha}} via \code{priors} argument in such a way:
#'
#'
#'\code{priors =list(mu_0 = c(1,1), S2 = ..., S3 = ..., alpha = ...)}
#'
#' with the constraint for \eqn{S2} and \eqn{S3} to be positive definite,
#' and \eqn{\bm{\alpha}} a vector of dimension \eqn{k} with nonnegative elements.
#'
#'When \code{software="rstan"}, the prior specification is:
#'
#'\deqn{ \bm{\mu}_j  \sim \mathcal{N}_D(\bm{\mu}_0, LD*L^{T})}
#'\deqn{L \sim \mbox{LKJ}(\epsilon)}
#'\deqn{D^*_j \sim \mbox{HalfCauchy}(0, \sigma_d).}
#'
#'The covariance matrix is expressed in terms of the LDL decomposition as \eqn{LD*L^{T}},
#'a variant of the classical Cholesky decomposition, where \eqn{L} is a \eqn{D \times D}
#'lower unit triangular matrix and \eqn{D*} is a \eqn{D \times D} diagonal matrix.
#'The Cholesky correlation factor \eqn{L} is assigned a LKJ prior with \eqn{\epsilon} degrees of freedom,  which,
#'combined with priors on the standard deviations of each component, induces a prior on the covariance matrix;
#'as \eqn{\epsilon \rightarrow \infty} the magnitude of correlations between components decreases,
#'whereas \eqn{\epsilon=1} leads to a uniform prior distribution for \eqn{L}.
#'By default, the hyperparameters are \eqn{\bm{\mu}_0=\bm{0}}, \eqn{\sigma_d=2.5, \epsilon=1}.
#'The user may propose some different values with the argument:
#'
#'
#' \code{priors=list(mu_0=c(1,2), sigma_d = 4, epsilon =2)}
#'
#'
#' If \code{software="rjags"} the function performs JAGS sampling using the \code{bayesmix} package
#' for univariate Gaussian mixtures, and the \code{runjags}
#' package for multivariate Gaussian mixtures. If \code{software="rstan"} the function performs
#' Hamiltonian Monte Carlo (HMC) sampling via the \code{rstan} package (see the vignette and the Stan project
#' for any help).
#'
#' After MCMC sampling, this function
#' clusters the units in \code{k} groups,
#' calls the \code{piv_sel()} function and yields the
#' pivots obtained from one among four different
#' methods (the user may specify one among them via \code{piv.criterion}
#' argument):
#'  \code{"maxsumint"}, \code{"minsumnoint"}, \code{"maxsumdiff"}
#'  and \code{"MUS"} (available only if \code{k <= 4})
#' (see the vignette for thorough details). Due to computational reasons
#' clarified in the Details section of the function \code{piv_rel}, the
#' length of the MCMC chains will be minor or equal than the input
#' argument \code{nMC}; this length, corresponding to the value
#' \code{true.iter} returned by the procedure, is the number of
#' MCMC iterations for which
#' the number of JAGS/Stan groups exactly coincides with the prespecified
#' number of groups \code{k}.
#' @return The function gives the MCMC output, the clustering
#' solutions and the pivotal indexes. Here there is a complete list of outputs.
#'
#' \item{\code{true.iter}}{ The number of MCMC iterations for which
#' the number of JAGS/Stan groups exactly coincides with the prespecified
#' number of groups \code{k}.}
#' \item{\code{Mu}}{An estimate of the groups' means.}
#' \item{\code{groupPost}}{ \eqn{true.iter \times N} matrix
#' with values from \code{1:k} indicating the post-processed group allocation
#' vector.}
#' \item{\code{mcmc_mean}}{  If \code{y} is a vector, a \eqn{true.iter \times k}
#' matrix with the post-processed MCMC chains for the mean parameters; if
#' \code{y} is a matrix, a \eqn{true.iter \times D \times k} array with
#' the post-processed MCMC chains for the mean parameters.}
#' \item{\code{mcmc_sd}}{  If \code{y} is a vector, a \eqn{true.iter \times k}
#' matrix with the post-processed MCMC chains for the sd parameters; if
#' \code{y} is a matrix, a \eqn{true.iter \times D} array with
#' the post-processed MCMC chains for the sd parameters.}
#' \item{\code{mcmc_weight}}{A \eqn{true.iter \times k}
#' matrix with the post-processed MCMC chains for the weights parameters.}
#'\item{\code{mcmc_mean_raw}}{ If \code{y} is a vector, a \eqn{(nMC-burn) \times k} matrix
#' with the raw MCMC chains for the mean parameters as given by JAGS; if
#' \code{y} is a matrix, a \eqn{(nMC-burn) \times D \times k} array with the raw MCMC chains
#' for the mean parameters as given by JAGS/Stan.}
#' \item{\code{mcmc_sd_raw}}{ If \code{y} is a vector, a \eqn{(nMC-burn) \times k} matrix
#' with the raw MCMC chains for the sd parameters as given by JAGS/Stan; if
#' \code{y} is a matrix, a \eqn{(nMC-burn) \times D} array with the raw MCMC chains
#' for the sd parameters as given by JAGS/Stan.}
#' \item{\code{mcmc_weight_raw}}{A \eqn{(nMC-burn) \times k} matrix
#' with the raw MCMC chains for the weights parameters as given by JAGS/Stan.}
#' \item{\code{C}}{The \eqn{N \times N} co-association matrix constructed from the MCMC sample.}
#' \item{\code{grr}}{The vector of cluster membership returned by
#' \code{"diana"} or \code{"hclust"}.}
#' \item{\code{pivots}}{The vector of indices of pivotal units identified by the selected pivotal criterion.}
#' \item{\code{model}}{The JAGS/Stan model code. Apply the \code{"cat"} function for a nice visualization of the code.}
#' \item{\code{stanfit}}{An object of S4 class \code{stanfit} for the fitted model (only if \code{software="rstan"}).}
#'
#' @author Leonardo Egidi \email{legidi@units.it}
#' @references Egidi, L., Pappadà, R., Pauli, F. and Torelli, N. (2018). Relabelling in Bayesian Mixture
#'Models by Pivotal Units. Statistics and Computing, 28(4), 957-969.
#' @examples
#'
#' ### Bivariate simulation
#'
#'\dontrun{
#' N   <- 200
#' k   <- 4
#' D   <- 2
#' nMC <- 1000
#' M1  <- c(-.5,8)
#' M2  <- c(25.5,.1)
#' M3  <- c(49.5,8)
#' M4  <- c(63.0,.1)
#' Mu  <- rbind(M1,M2,M3,M4)
#' Sigma.p1 <- diag(D)
#' Sigma.p2 <- 20*diag(D)
#' W <- c(0.2,0.8)
#' sim <- piv_sim(N = N, k = k, Mu = Mu,
#'                Sigma.p1 = Sigma.p1,
#'                Sigma.p2 = Sigma.p2, W = W)
#'
#' ## rjags (default)
#' res <- piv_MCMC(y = sim$y, k =k, nMC = nMC)
#'
#' ## rstan
#' res_stan <- piv_MCMC(y = sim$y, k =k, nMC = nMC,
#'                      software ="rstan")
#'
#' # changing priors
#' res2 <- piv_MCMC(y = sim$y,
#'                  priors = list (
#'                  mu_0=c(1,1),
#'                  S2 = matrix(c(0.002,0,0, 0.1),2,2, byrow=TRUE),
#'                  S3 = matrix(c(0.1,0,0,0.1), 2,2, byrow =TRUE)),
#'                  k = k, nMC = nMC)
#'}
#'
#'
#' ### Fishery data (bayesmix package)
#'
#'\dontrun{
#' library(bayesmix)
#' data(fish)
#' y <- fish[,1]
#' k <- 5
#' nMC <- 5000
#' res <- piv_MCMC(y = y, k = k, nMC = nMC)
#'
#' # changing priors
#' res2   <- piv_MCMC(y = y,
#'                    priors = list(kind = "condconjugate",
#'                    parameter = "priorsRaftery",
#'                    hierarchical = "tau"),  k =k, nMC = nMC)
#'}
#' @export






piv_MCMC <- function(y,
                     k,
                     nMC,
                     priors,
                     piv.criterion = c("MUS", "maxsumint", "minsumnoint", "maxsumdiff"),
                     clustering = c("diana", "hclust"),
                     software = c("rjags", "rstan"),
                     burn =0.5*nMC,
                     chains = 4,
                     cores = 1,
                     sparsity = TRUE){

  #### checks

  # piv.criterion
  if (missing(piv.criterion)){
    piv.criterion <- "maxsumint"
  }
  list_crit <- c("MUS", "maxsumint", "minsumnoint", "maxsumdiff")
  piv.criterion <- match.arg(piv.criterion, list_crit)

  # clustering

  list_clust <- c("diana", "hclust")
  clustering <- match.arg(clustering, list_clust)

  # software

  list_soft <- c("rjags", "rstan")
  software <- match.arg(software, list_soft)

  # burn-in

  if (burn > nMC){
    stop("Please, 'burn' argument has to be minor than
         the number of MCMC iterations!", sep="")
  }


  ###

  # Conditions about data dimension----------------
  if (missing(software)){
    software="rjags"
  }

  if (is.vector(y)){
    N <- length(y)
    # Initial values
    mu_inits<- c()
    clust_inits <- kmeans(y, k, nstart = 10)$cluster
    for (j in 1:k){
      mu_inits[j]<-mean(y[clust_inits==j])
    }
    if (software=="rjags"){
      priors_input <- list(kind = "independence",
                           parameter = "priorsFish",
                           hierarchical = "tau")
      if (missing(priors)){
        # b0 = 0; B0inv =0.1; nu0Half =5;
        # g0Half = 1e-17; g0G0Half = 1e16;
        # e = rep(1,k); S0 =2
        # priors =  list( "b0" , "B0inv" , "nu0Half",
        #                 "g0Half", "g0G0Half", "e", "S0")

        mod.mist.univ <- BMMmodel(y, k = k,
                                  initialValues = list(S0 = 2),
                                  priors = priors_input)
      }else{
        if (is.null(priors$mu_0)){
          b0 <- median(as.matrix(y))
        }else{
          b0 <- priors$mu_0
        }

        if (is.null(priors$B0inv)){
          B0inv <- 1/priorsFish(y)$B0
        }else{
          B0inv <- priors$B0inv
        }

        if (is.null(priors$nu_0)){
          nu0Half <- priorsFish(y)$nu0/2
        }else{
          nu0Half <- priors$nu_0/2
        }

        if (is.null(priors$g_0)){
          g0Half <- 0.5*10^-16
        }else{
          g0Half <- priors$g_0/2
        }

        if (is.null(priors$G_0)){
          g0G0Half <- 0.5*10^-16
        }else{
          g0G0Half <- priors$G_0/2
        }

        if (is.null(priors$alpha)){
          e <- rep(0.001,k)
        }else{
          e <- priors$alpha
        }

        if (is.null(priors$S0)){
          S0 <- priorsFish(y)$S0
        }else{
          S0 <- priors$S0
        }


        nu0S0Half = nu0Half*S0
        mod.mist.univ <- BMMmodel(y, k = k,
                                  initialValues = list(S0 = 2),
                                  priors = priors_input)
        #prior's values redefinition
        mod.mist.univ$data$b0 <- b0
        mod.mist.univ$data$B0inv <- B0inv
        mod.mist.univ$data$nu0Half <- nu0Half
        mod.mist.univ$data$g0Half <- g0Half
        mod.mist.univ$data$g0G0Half <- g0G0Half
        mod.mist.univ$data$e <- e

      }

      # JAGS code------------------------

      # Data
      # Model

      control <- JAGScontrol(variables = c("mu", "tau", "eta", "S"),
                             burn.in = burn, n.iter = nMC, seed = 10)
      ogg.jags <- JAGSrun(y, model = mod.mist.univ, control = control)
      printed <- cat(print(ogg.jags))
      # Parameters' initialization

      J <- 3
      mcmc.pars <- array(data = NA, dim = c(nMC-length(1:burn), k, J))
      mcmc.pars[ , , 1] <- ogg.jags$results[-(1:burn), (N+k+1):(N+2*k)]
      mcmc.pars[ , , 2] <- ogg.jags$results[-(1:burn), (N+2*k+1):(N+3*k)]
      mcmc.pars[ , , 3] <- ogg.jags$results[-(1:burn), (N+1):(N+k)]

      mu_pre_switch_compl <-  mcmc.pars[ , , 1]
      tau_pre_switch_compl <-  1/mcmc.pars[ , , 2]
      prob.st_pre_switch_compl <-  mcmc.pars[ , , 3]

      mu <- mcmc.pars[,,1]
      tau <- 1/mcmc.pars[,,2]
      prob.st <- mcmc.pars[,,3]
      group <-  ogg.jags$results[-(1:burn), 1:N] #gruppi
      group_for_nclusters <- group
      FreqGruppiJags <- table(group)
      numeffettivogruppi <- apply(group,1,FUN = function(x) length(unique(x)))

      if (sum(numeffettivogruppi==k)==0){
        print("MCMC has not never been able to identify the required number of groups and the process has been interrupted")
        return(list(nclusters = numeffettivogruppi))

        #return(1)
      }

      ##saved in the output
      ris_prel <- ogg.jags$results[-(1:burn),]
      ris <- ris_prel[numeffettivogruppi==k,]
      true.iter <- nrow(ris)
      group <- ris[,1:N]
      mu <- mu[numeffettivogruppi==k,]
      tau <- tau[numeffettivogruppi==k,]
      prob.st <- prob.st[numeffettivogruppi==k,]
      model_code <- mod.mist.univ$bugs
      ## just another way to save them
      mcmc_mean_raw <- mcmc.pars[,,1]
      mcmc_sd_raw <- mcmc.pars[,,2]
      mcmc_weight_raw <- mcmc.pars[,,3]

    }else if (software=="rstan"){
      if(missing(priors)){
        mu_0 <- 0
        B0inv <- 0.1
        mu_sigma <- 0
        tau_sigma <- 2
        a_sp <- 1
        b_sp <- 200
      }else{
        if (is.null(priors$mu_0)){
          mu_0 <- 0
        }else{
          mu_0 <- priors$mu_0
        }
        if (is.null(priors$B0inv)){
          B0inv <- 0.1
        }else{
          B0inv <- priors$B0inv
        }
        if (is.null(priors$mu_sigma)){
          mu_sigma <- 0
        }else{
          mu_sigma <- priors$mu_sigma
        }
        if (is.null(priors$tau_sigma)){
          tau_sigma <- 2
        }else{
          tau_sigma <- priors$tau_sigma
        }
        if (is.null(priors$a_sp)){
          a_sp <- 1
        }else{
          a_sp <- priors$a_sp
        }
        if (is.null(priors$b_sp)){
          b_sp <- 200
        }else{
          b_sp <- priors$b_sp
        }
      }

      data = list(N=N, y=y, k=k,
                  mu_0=mu_0, B0inv=B0inv,
                  mu_sigma=mu_sigma, tau_sigma=tau_sigma,
                  a=a_sp, b=b_sp)
      # sparsity
      if (sparsity == FALSE){

      mix_univ <-"
        data {
          int<lower=1> k;          // number of mixture components
          int<lower=1> N;          // number of data points
          real y[N];               // observations
          real mu_0;               // mean hyperparameter
          real<lower=0> B0inv;     // mean hyperprecision
          real mu_sigma;           // sigma hypermean
          real<lower=0> tau_sigma; // sigma hyper sd
          }
        parameters {
          simplex[k] eta;             // mixing proportions
          ordered[k] mu;              // locations of mixture components
          vector<lower=0>[k] sigma;   // scales of mixture components
          }
      transformed parameters{
          vector[k] log_eta = log(eta);  // cache log calculation
          vector[k] pz[N];
          simplex[k] exp_pz[N];
              for (n in 1:N){
                  pz[n] =   normal_lpdf(y[n]|mu, sigma)+
                            log_eta-
                            log_sum_exp(normal_lpdf(y[n]|mu, sigma)+
                            log_eta);
                  exp_pz[n] = exp(pz[n]);
                            }
          }
      model {
        sigma ~ lognormal(mu_sigma, tau_sigma);
        mu ~ normal(mu_0, 1/B0inv);
            for (n in 1:N) {
              vector[k] lps = log_eta;
                for (j in 1:k){
                    lps[j] += normal_lpdf(y[n] | mu[j], sigma[j]);
                    target+=pz[n,j];
                    }
              target += log_sum_exp(lps);
                  }
          }
     generated quantities{
        int<lower=1, upper=k> z[N];
          for (n in 1:N){
              z[n] = categorical_rng(exp_pz[n]);
            }
      }
      "
      }else{
      mix_univ <-"
      data {
        int<lower=1> k;          // number of mixture components
        int<lower=1> N;          // number of data points
        real y[N];               // observations
        real mu_0;               // mean hyperparameter
        real<lower=0> B0inv;     // mean hyperprecision
        real mu_sigma;           // sigma hypermean
        real<lower=0> tau_sigma; // sigma hyper sd
        real<lower=0> a;         // hyper-shape gamma e0
        real<lower=0> b;         // hyper-rate gamma e0
      }
      parameters {
        simplex[k] eta;             // mixing proportions
        ordered[k] mu;              // locations of mixture components
        vector<lower=0>[k] sigma;   // scales of mixture components
        real<lower=0> e0;
      }
      transformed parameters{
        vector[k] log_eta = log(eta);  // cache log calculation
        vector<lower=0>[k] alpha = rep_vector(e0, k);
        vector[k] pz[N];
        simplex[k] exp_pz[N];
        for (n in 1:N){
          pz[n] =   normal_lpdf(y[n]|mu, sigma)+
            log_eta-
            log_sum_exp(normal_lpdf(y[n]|mu, sigma)+
                          log_eta);
          exp_pz[n] = exp(pz[n]);
        }
      }
      model {
        sigma ~ lognormal(mu_sigma, tau_sigma);
        mu ~ normal(mu_0, 1/B0inv);
        eta ~ dirichlet(alpha);
        e0 ~ gamma(a, b);
        for (n in 1:N) {
          vector[k] lps = log_eta;
          for (j in 1:k){
            lps[j] += normal_lpdf(y[n] | mu[j], sigma[j]);
            target+=pz[n,j];
          }
          target += log_sum_exp(lps);
        }
      }
      generated quantities{
        int<lower=1, upper=k> z[N];
        for (n in 1:N){
          z[n] = categorical_rng(exp_pz[n]);
        }
      }
      "
      }

      fit_univ <-  stan(model_code = mix_univ,
                        data=data,
                        chains =chains,
                        iter =nMC)
      stanfit <- fit_univ
      if (sparsity==FALSE){
        printed <- cat(print(fit_univ, pars =c("mu", "eta", "sigma")))
      }else{
        printed <- cat(print(fit_univ, pars =c("mu", "eta", "sigma", "e0")))
      }
      sims_univ <- rstan::extract(fit_univ)

      J <- 3
      mcmc.pars <- array(data = NA, dim = c(dim(sims_univ$eta)[1], k, J))
      mcmc.pars[ , , 1] <- sims_univ$mu
      mcmc.pars[ , , 2] <- sims_univ$sigma
      mcmc.pars[ , , 3] <- sims_univ$eta

      mu_pre_switch_compl <-  mcmc.pars[ , , 1]
      tau_pre_switch_compl <-  mcmc.pars[ , , 2]
      prob.st_pre_switch_compl <-  mcmc.pars[ , , 3]

      mu <- mcmc.pars[,,1]
      tau <- mcmc.pars[,,2]
      prob.st <- mcmc.pars[,,3]
      group <-  sims_univ$z[, 1:N] #gruppi
      group_for_nclusters <- group
      FreqGruppiJags <- table(group)
      numeffettivogruppi <- apply(group,1,FUN = function(x) length(unique(x)))

      if (sum(numeffettivogruppi==k)==0){
        print("MCMC has not never been able to identify the required number of groups and the process has been interrupted")
        return(list(nclusters = numeffettivogruppi))

      }

      ##saved in the output
      ris_prel <- as.matrix(fit_univ)
      ris <- ris_prel[numeffettivogruppi==k,]
      group <- group[numeffettivogruppi==k,]
      mu <- mu[numeffettivogruppi==k,]
      tau <- tau[numeffettivogruppi==k,]
      prob.st <- prob.st[numeffettivogruppi==k,]
      true.iter <- nrow(ris)
      model_code <- mix_univ
      ## just another way to save them
      mcmc_mean_raw <- mcmc.pars[,,1]
      mcmc_sd_raw <- mcmc.pars[,,2]
      mcmc_weight_raw <- mcmc.pars[,,3]

    }


    ## resambling

    group.orig <- group
    verigruppi <- as.double(names(table(group)))

    cont <- 0
    for (verogruppo in verigruppi){
      cont <- cont+1
      # update pivot counter
      group.orig[group==verogruppo] <- cont
    }
    cont

    k.orig <- k
    if (cont>1){
      k <- cont
    }
    mu <- mu[,verigruppi]
    tau <- tau[,verigruppi]
    prob.st <- prob.st[,verigruppi]

    M <- nrow(group)
    group <- group*0
    mu_switch = tau_switch = prob.st_switch = array(rep(0, true.iter*k), dim=c(true.iter,k))
    z <- array(0,dim=c(N, k, true.iter))

    for (i in 1:true.iter){
      perm <- sample(1:k,k,replace=FALSE)
      for (j in 1:k){
        #post-processing
        group[i,group.orig[i,]==j] <- perm[j]
      }
      mu_switch[i,] <- mu[i,perm]
      tau_switch[i,] <- tau[i,perm]
      prob.st_switch[i,] <- prob.st[i,perm]
    }

    for (i in 1:true.iter){
      for (j in 1:N){
        z[j,group[i,j],i] <- 1
      }
    }

    mcmc_mean = mcmc_sd = mcmc_weight = array(0, dim=c(true.iter, k))
    mcmc_mean <- mu_switch
    mcmc_sd <- tau_switch
    mcmc_weight <- prob.st_switch



  }else if (is.matrix(y)){
    N <- dim(y)[1]
    D <- dim(y)[2]
    # Parameters' initialization
    clust_inits <- kmeans(y, k, nstart = 10)$cluster
    mu_inits <- matrix(0,k,D)
    for (j in 1:k){
      for (d in 1:D){
      mu_inits[j, d] <- mean(y[clust_inits==j,d])
    }}
    #Reorder mu_inits according to the x-coordinate
    mu_inits <-
      mu_inits[sort(mu_inits[,1], decreasing=FALSE, index.return=TRUE)$ix,]

    if (software=="rjags"){

      # JAGS code------------------------

      # Initial values
      if (missing(priors)){
        mu_0 <- rep(0, D)
        S2 <- diag(D)/100000
        S3 <- diag(D)/100000
        alpha <- rep(1,k)
        a_sp <-1
        b_sp <- 200
      }else{
        if (is.null(priors$mu_0)){
          mu_0 <- rep(0, D)
        }else{
          mu_0 <- priors$mu_0
        }
        if (is.null(priors$S2)){
          S2 <- diag(D)/100000
        }else{
          S2 <- priors$S2
        }
        if (is.null(priors$S3)){
          S3 <- diag(D)/100000
        }else{
          S3 <- priors$S3
        }
        if (is.null(priors$alpha)){
          alpha <- rep(1, k)
        }else{
          alpha <- priors$alpha
        }
        if (is.null(priors$a_sp)){
          a_sp <- 1
        }else{
          a_sp <- priors$a_sp
        }
        if (is.null(priors$b_sp)){
          b_sp <- 200
        }else{
          b_sp <- priors$b_sp
        }
      }

      # Model
      if (sparsity ==FALSE){

        # Data
        dati.biv <- list(y = y, N = N, k = k, D = D,
                         S2= S2, S3= S3, mu_0=mu_0,
                         alpha = alpha)
      mod.mist.biv<-"model{
    # Likelihood:

    for (i in 1:N){
      yprev[i,1:D]<-y[i,1:D]
      y[i,1:D] ~ dmnorm(mu[clust[i],],tau)
      clust[i] ~ dcat(eta[1:k] )
    }

    # Prior:

    for (g in 1:k) {
      mu[g,1:D] ~ dmnorm(mu_0[],S2[,])}
      tau[1:D,1:D] ~ dwish(S3[,],D+1)
      Sigma[1:D,1:D] <- inverse(tau[,])
      eta[1:k] ~ ddirch(alpha)
  }"
      }else{

        # Data
        dati.biv <- list(y = y, N = N, k = k, D = D,
                         S2= S2, S3= S3, mu_0=mu_0,
                         a = a_sp, b=b_sp)
        mod.mist.biv<-"model{
      # Likelihood:

      for (i in 1:N){
        yprev[i,1:D]<-y[i,1:D]
        y[i,1:D] ~ dmnorm(mu[clust[i],],tau)
        clust[i] ~ dcat(eta[1:k] )
      }

      # Prior:

      for (g in 1:k) {
        mu[g,1:D] ~ dmnorm(mu_0[],S2[,])
        alpha[g] <- e0}
        tau[1:D,1:D] ~ dwish(S3[,],D+1)
        Sigma[1:D,1:D] <- inverse(tau[,])
        eta[1:k] ~ ddirch(alpha)
        e0 ~ dgamma(a, b)
    }"

      }

      init1.biv <- list()
      for (s in 1:chains){
      init1.biv[[s]] <- dump.format(list(mu=mu_inits,
                                    tau= 15*diag(D),
                                    eta=rep(1/k,k), clust=clust_inits))
      }
      moni.biv <- c("clust","mu","tau","eta", "e0")

      mod   <- mod.mist.biv
      dati  <- dati.biv
      init1 <- init1.biv
      moni  <- moni.biv

      # Jags execution
      ogg.jags <- run.jags(model=mod, data=dati, monitor=moni,
                           inits=init1, n.chains=chains,plots=FALSE, thin=1,
                           sample=nMC, burnin=burn)
      printed <- print(add.summary(ogg.jags, vars= c("mu",
                      "tau",
                      "eta")))
      # Extraction
      ris <- ogg.jags$mcmc[[1]]

      # Post- process of the chains----------------------
      group <- ris[-(1:burn),grep("clust[",colnames(ris),fixed=TRUE)]
      group_for_nclusters <- group
      # only the variances
      tau <- sqrt( (1/ris[-(1:burn),grep("tau[",colnames(ris),fixed=TRUE)])[,c(1,4)])
      prob.st <- ris[-(1:burn),grep("eta[",colnames(ris),fixed=TRUE)]
      M <- nrow(group)
      H <- list()

      mu_pre_switch_compl <- array(rep(0, M*D*k), dim=c(M,D,k))
      for (i in 1:k){
        H[[i]] <- ris[-(1:burn),
                      grep(paste("mu[",i, sep=""),
                           colnames(ris),fixed=TRUE)]
                         #[,c(i,i+k)]
      }
      for (i in 1:k){
        mu_pre_switch_compl[,,i] <- as.matrix(H[[i]])
      }
      # Discard iterations
      numeffettivogruppi <- apply(group,1,FUN = function(x) length(unique(x)))
      ris <- ris[numeffettivogruppi==k,]
      true.iter <- nrow(ris)

      if (sum(numeffettivogruppi==k)==0){
        print("MCMC has not never been able to identify the required number of groups and the process has been interrupted")
        return(list(nclusters = numeffettivogruppi))
      }else{
        L<-list()
        mu_pre_switch <- array(rep(0, true.iter*D*k), dim=c(true.iter,D,k))
        for (i in 1:k){
          L[[i]] <- ris[,grep(paste("mu[", i, sep=""),
                                    colnames(ris),fixed=TRUE)]
        }
        for (i in 1:k){
          mu_pre_switch[,,i] <- as.matrix(L[[i]])
        }
      }

      group <- ris[,grep("clust[",colnames(ris),fixed=TRUE)]
      FreqGruppiJags <- table(group)
      model_code <- mod.mist.biv

      mcmc_mean_raw = mu_pre_switch_compl
      mcmc_weight_raw = prob.st
      mcmc_sd_raw = tau

      tau <- sqrt( (1/ris[,grep("tau[",colnames(ris),fixed=TRUE)])[,c(1,4)])
      prob.st <- ris[,grep("eta[",colnames(ris),fixed=TRUE)]
      mu <- mu_pre_switch


    }else if(software=="rstan"){
      if (missing(priors)){
        mu_0 <- rep(0, D)
        epsilon <- 1
        sigma_d <- 2.5
        a_sp <- 1
        b_sp <- 200
      }else{
        if (is.null(priors$mu_0)){
          mu_0 <- rep(0, D)
        }else{
          mu_0 <- priors$mu_0
        }
        if (is.null(priors$epsilon)){
          epsilon <- 1
        }else{
          epsilon  <- priors$epsilon
        }
        if (is.null(priors$sigma_d)){
          sigma_d <- 2.5
        }else{
          sigma_d <- priors$sigma_d
        }
        if (is.null(priors$a_sp)){
          a_sp <- 1
        }else{
          a_sp <- priors$a_sp
        }
        if (is.null(priors$b_sp)){
          b_sp <- 200
        }else{
          b_sp <- priors$b_sp
        }
      }
      data =list(N=N, k=k, y=y, D=D, mu_0=mu_0,
                 epsilon = epsilon, sigma_d = sigma_d,
                 a = a_sp, b=b_sp)
      # sparsity
      if (sparsity=="FALSE"){
      mix_biv <- "
        data {
          int<lower=1> k;          // number of mixture components
          int<lower=1> N;          // number of data points
          int D;                   // data dimension
          matrix[N,D] y;           // observations matrix
          vector[D] mu_0;
          real<lower=0> epsilon;
          real<lower=0> sigma_d;
        }
        parameters {
          simplex[k] eta;         // mixing proportions
          vector[D] mu[k];        // locations of mixture components
          cholesky_factor_corr[D] L_Omega;   // scales of mixture components
          vector<lower=0>[D] L_sigma;
          cholesky_factor_corr[D] L_tau_Omega;   // scales of mixture components
          vector<lower=0>[D] L_tau;
          }
        transformed parameters{
          vector[k] log_eta = log(eta);  // cache log calculation
          vector[k] pz[N];
          simplex[k] exp_pz[N];
          matrix[D,D] L_Sigma=diag_pre_multiply(L_sigma, L_Omega);
          matrix[D,D] L_Tau=diag_pre_multiply(L_tau, L_tau_Omega);


            for (n in 1:N){
                pz[n]=   multi_normal_cholesky_lpdf(y[n]|mu, L_Sigma)+
                         log_eta-
                         log_sum_exp(multi_normal_cholesky_lpdf(y[n]|
                                                     mu, L_Sigma)+
                         log_eta);
                exp_pz[n] = exp(pz[n]);
              }
          }
        model{
          L_Omega ~ lkj_corr_cholesky(epsilon);
          L_sigma ~ cauchy(0, sigma_d);
          mu ~ multi_normal_cholesky(mu_0, L_Tau);
            for (n in 1:N) {
              vector[k] lps = log_eta;
                for (j in 1:k){
                    lps[j] += multi_normal_cholesky_lpdf(y[n] |
                                                   mu[j], L_Sigma);
                    target+=pz[n,j];
                }
              target += log_sum_exp(lps);
          }
          }
        generated quantities{
          int<lower=1, upper=k> z[N];
              for (n in 1:N){
                  z[n] = categorical_rng(exp_pz[n]);
                }
          }
          "
      }else{
        mix_biv <- "
        data {
          int<lower=1> k;          // number of mixture components
          int<lower=1> N;          // number of data points
          int D;                   // data dimension
          matrix[N,D] y;           // observations matrix
          vector[D] mu_0;
          real<lower=0> epsilon;
          real<lower=0> sigma_d;
          real<lower=0> a;         // hyper-shape gamma e0
          real<lower=0> b;         // hyper-rate gamma e0
        }
        parameters {
          simplex[k] eta;         // mixing proportions
          vector[D] mu[k];        // locations of mixture components
          cholesky_factor_corr[D] L_Omega;   // scales of mixture components
          vector<lower=0>[D] L_sigma;
          cholesky_factor_corr[D] L_tau_Omega;   // scales of mixture components
          vector<lower=0>[D] L_tau;
          real<lower=0> e0;       // dirichlet concentration
          }
        transformed parameters{
          vector[k] log_eta = log(eta);  // cache log calculation
          vector<lower=0>[k] alpha = rep_vector(e0, k);
          vector[k] pz[N];
          simplex[k] exp_pz[N];
          matrix[D,D] L_Sigma=diag_pre_multiply(L_sigma, L_Omega);
          matrix[D,D] L_Tau=diag_pre_multiply(L_tau, L_tau_Omega);


            for (n in 1:N){
                pz[n]=   multi_normal_cholesky_lpdf(y[n]|mu, L_Sigma)+
                         log_eta-
                         log_sum_exp(multi_normal_cholesky_lpdf(y[n]|
                                                     mu, L_Sigma)+
                         log_eta);
                exp_pz[n] = exp(pz[n]);
              }
          }
        model{
          L_Omega ~ lkj_corr_cholesky(epsilon);
          L_sigma ~ cauchy(0, sigma_d);
          mu ~ multi_normal_cholesky(mu_0, L_Tau);
          eta ~ dirichlet(alpha);
          e0 ~ gamma(a, b);
            for (n in 1:N) {
              vector[k] lps = log_eta;
                for (j in 1:k){
                    lps[j] += multi_normal_cholesky_lpdf(y[n] |
                                                   mu[j], L_Sigma);
                    target+=pz[n,j];
                }
              target += log_sum_exp(lps);
          }
          }
        generated quantities{
          int<lower=1, upper=k> z[N];
              for (n in 1:N){
                  z[n] = categorical_rng(exp_pz[n]);
                }
          }
          "


      }
      fit_biv <-  stan(model_code = mix_biv,
                       data=data,
                       chains =chains,
                       iter =nMC)
      stanfit <- fit_biv
      if (sparsity == FALSE){
        printed <- cat(print(fit_biv, pars=c("mu", "eta", "L_Sigma")))
      }else{
        printed <- cat(print(fit_biv, pars=c("mu", "eta", "L_Sigma", "e0")))
      }
      sims_biv <- rstan::extract(fit_biv)

      # Extraction
      ris <- as.matrix(sims_biv)

      # Post- process of the chains----------------------
      group <- sims_biv$z
      group_for_nclusters <- group
      tau <- sims_biv$L_sigma
      prob.st <- sims_biv$eta
      M <- nrow(group)

      mu_pre_switch_compl <- array(rep(0, M*D*k), dim=c(M,D,k))
      for (i in 1:M)
        mu_pre_switch_compl[i,,] <- t(sims_biv$mu[i,,])
      # Discard iterations
      numeffettivogruppi <- apply(group,1,FUN = function(x) length(unique(x)))
      sm <- sims_biv$mu[numeffettivogruppi==k,,]
      true.iter <- dim(sm)[1]

      if (sum(numeffettivogruppi==k)==0){
        print("MCMC has not never been able to identify the required number of groups and the process has been interrupted")
        return(list(nclusters = numeffettivogruppi))
      }else{

        mu_pre_switch <- array(rep(0, true.iter*D*k), dim=c(true.iter,D,k))
        for (i in 1:true.iter)
          mu_pre_switch[i,,] <- t(sm[i,,])
      }

      mcmc_mean_raw = mu_pre_switch_compl
      mcmc_weight_raw = prob.st
      mcmc_sd_raw = tau

      group <- sims_biv$z[numeffettivogruppi==k,]
      mu <- mu_pre_switch
      tau <- sims_biv$L_sigma[numeffettivogruppi==k, ]
      prob.st <- sims_biv$eta[numeffettivogruppi==k,]
      FreqGruppiJags <- table(group)

      model_code <- mix_biv
    }

    group.orig <- group
    verigruppi <- as.double(names(table(group)))
    prob.st <- prob.st[,verigruppi]
    mu <- mu[,,verigruppi]

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
    mu_switch  <- array(rep(0, true.iter*D*k), dim=c(true.iter,D,k))
    prob.st_switch <-  array(0, dim=c(true.iter,k))
    group <- group*0
    z <- array(0,dim=c(N, k, true.iter))

    for (i in 1:true.iter){
      perm <- sample(1:k,k,replace=FALSE)
      for (j in 1:k){
        #post-processing
        group[i,group.orig[i,]==j] <- perm[j]
      }
      mu_switch[i,,] <- mu[i,,perm]
      prob.st_switch[i,] <- prob.st[i,perm]
    }

    for (i in 1:true.iter){
      for (j in 1:N){
        z[j,group[i,j],i] <- 1
      }
    }

    mcmc_mean <- mu_switch
    mcmc_sd <- tau
    mcmc_weight <- prob.st_switch


  }

  FreqGruppiJagsPERM <- table(group)
  Freq <- cbind(FreqGruppiJags,FreqGruppiJagsPERM)
  colnames(Freq) <- c("JAGS raw groups", "JAGS post-processed groups")



  # Similarity matrix based on MCMC sampling------------------------
  nz <- dim(z)[1]
  M <- dim(z)[3]
  C <- matrix(NA,nz,nz)
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

  piv.criterion.choices <- c("maxsumint", "minsumnoint",
                             "maxsumdiff")

  if (missing(piv.criterion)){
    piv.criterion <- "maxsumint"
  }

  if (piv.criterion=="maxsumint"||
      piv.criterion=="minsumnoint"||
      piv.criterion=="maxsumdiff" ){

    piv.index <- (1:3)[piv.criterion.choices==piv.criterion]
    piv.index.pivotal <- c(1,2,3)
    available_met <- 3
    x <- c(1:available_met)
    clust  <-  piv_sel(C=C, clusters=as.vector(grr))
    pivots <- clust$pivots[,piv.index.pivotal[piv.index]]
  }else if(piv.criterion=="MUS"){
    if (k <=4 & sum(C==0)!=0){

      mus_res    <- MUS(C, grr)
      pivots     <-  mus_res$pivots

    }else{

      print("maxsumint criterion instead of MUS has been adopted due to
          computational efficiency")
      clust  <-  piv_sel(C=C,  clusters=as.vector(grr))
      pivots <- clust$pivots[,3]
    }
  }

  if (software == "rjags"){
    stanfit = NULL
  }



  return(list( true.iter = true.iter,
               Mu = mu_inits,
               groupPost=group,
               mcmc_mean = mcmc_mean,
               mcmc_sd = mcmc_sd,
               mcmc_weight = mcmc_weight,
               mcmc_mean_raw = mcmc_mean_raw,
               mcmc_sd_raw = mcmc_sd_raw,
               mcmc_weight_raw = mcmc_weight_raw,
               C=C,
               grr=grr,
               pivots = pivots,
               model = model_code,
               k = k,
               stanfit = stanfit,
               nclusters = numeffettivogruppi))
}
