#' Selection of the number of mixture components k via
#' leave-one-out (loo) cross-validation.
#'
#'
#'
#'
#'
#'
#'
#'
#' @export

piv_k <- function(y, priors,  ic =c("dic", "loo"),
                  nMC = 500, chains = 4, k_max = 8, burn = 0.5*nMC){
  if (is.vector(y)){
    N <- length(y)

    if (ic=="dic"){

    k_vec <- c(1:k_max)
    DIC <- c()

    mod.mist.univ  <- "model {
                      # Likelihood:
                        for( i in 1 : N ) {
                          y[i] ~ dnorm( mu[i] , tau[i] )
                          mu[i] <- muOfClust[ clust[i] ]
                          tau[i] <- tauOfClust[ clust[i]]
                          clust[i] ~ dcat( pClust[1:k] )
                                }
                       # Prior:
                          for ( g in 1: k ) {
                            tauOfClust[g] ~ dgamma( 0.01 , 0.01 )
                            muOfClust[g] ~ dnorm( mu_0 , 1.0E-10 )
                              }
                          pClust[1:k] ~ ddirch( alpha )
                              }
                              "


    for (h in 2:k_max){
      if (missing(priors)){
        mu_0 <- 0
        alpha <- rep(1,k_vec[h])
      }
      dati.univ <- list(y = y, N = N, k = k_vec[h],
                        mu_0= mu_0,
                       alpha = alpha)
      # Initial values
      mu_inits<- c()
      clust_inits <- kmeans(y, k_vec[h])$cluster
      for (j in 1:k_vec[h]){
        mu_inits[j]<-mean(y[clust_inits==j])
      }
      init1.univ <- dump.format(list(muOfClust=mu_inits,
                                     tauOfClust= runif(k_vec[h], 0,10),
                                     pClust=rep(1/k_vec[h],k_vec[h]),
                                     clust=clust_inits))
      init2.univ <- dump.format(list(muOfClust=mu_inits,
                                     tauOfClust= runif(k_vec[h], 0,10),
                                     pClust=rep(1/k_vec[h],k_vec[h]),
                                     clust=clust_inits))
      init3.univ <- dump.format(list(muOfClust=mu_inits,
                                     tauOfClust= runif(k_vec[h], 0,10),
                                     pClust=rep(1/k_vec[h],k_vec[h]),
                                     clust=clust_inits))

      moni.univ <- c("clust","muOfClust","tauOfClust","pClust", "dic")

      mod   <- mod.mist.univ
      dati  <- dati.univ
      init1 <- list(init1.univ, init2.univ, init3.univ)
      moni  <- moni.univ

      # Jags execution
      ogg.jags <- run.jags(model=mod, data=dati, monitor=moni,
                           inits=init1, n.chains=chains,plots=FALSE, thin=1,
                           sample=nMC, burnin=burn)
      ris <- ogg.jags$mcmc[[1]]
      DIC[h-1] <- as.double(ogg.jags$dic$dic)
    }
    k_est <- which.min(DIC)+1
    return(list(dic=DIC, k_est =k_est))
    }else if (ic == "loo"){
      if(missing(priors)){
        mu_0 <- 0
        B0inv <- 0.1
        mu_phi <- 0
        sigma_phi <- 2
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
        if (is.null(priors$mu_phi)){
          mu_phi <- 0
        }else{
          mu_phi <- priors$mu_phi
        }
        if (is.null(priors$sigma_phi)){
          sigma_phi <- 2
        }else{
          sigma_phi <- priors$sigma_phi
        }
      }


      mix_univ <-"
        data {
          int<lower=1> k;          // number of mixture components
          int<lower=1> N;          // number of data points
          real y[N];               // observations
          real mu_0;               // mean hyperparameter
          real<lower=0> B0inv;     // mean hyperprecision
          real mu_phi;             // sigma hypermean
          real<lower=0> sigma_phi; // sigma hyper sd
          }
        parameters {
          simplex[k] theta;        // mixing proportions
          ordered[k] mu;              // locations of mixture components
          vector<lower=0>[k] sigma;   // scales of mixture components
          }
      transformed parameters{
          vector[k] log_theta = log(theta);  // cache log calculation
          vector[k] pz[N];
          simplex[k] exp_pz[N];
              for (n in 1:N){
                  pz[n] =   normal_lpdf(y[n]|mu, sigma)+
                            log_theta-
                            log_sum_exp(normal_lpdf(y[n]|mu, sigma)+
                            log_theta);
                  exp_pz[n] = exp(pz[n]);
                            }
          }
      model {

        sigma ~ lognormal(mu_phi, sigma_phi);
        mu ~ normal(mu_0, 1/B0inv);
            for (n in 1:N) {
              vector[k] lps = log_theta;
                for (j in 1:k){
                    lps[j] += normal_lpdf(y[n] | mu[j], sigma[j]);
                    target+=pz[n,j];
                    }
              target += log_sum_exp(lps);
                  }
          }
     generated quantities{
        matrix[N,k] ll;
        vector[N] log_lik;
        int<lower=1, upper=k> z[N];
          for (n in 1:N){
              z[n] = categorical_rng(exp_pz[n]);
          }

        for (n in 1:N){
          for (j in 1:k){
            ll[n,j] = normal_lpdf(y[n] | mu[j], sigma[j])+log_theta[j];
          }
          log_lik[n] = log_sum_exp(ll[n,]);
        }
      }
      "

      k_vec <- c(1:k_max)
      LOOIC <- c()

      for (h in 2:k_max){
        data <- list(N=N, y=y, k=k_vec[h],
                    mu_0=mu_0, B0inv=B0inv,
                    mu_phi=mu_phi, sigma_phi=sigma_phi)
        fit_univ <-  stan(model_code = mix_univ,
                        data=data,
                        chains =chains,
                        iter =nMC)
        log_lik <- extract_log_lik(fit_univ)
        loo <- loo(log_lik)
        LOOIC[h-1] <- loo$estimates[3,1]
      }
      k_est <- which.min(LOOIC)+1
     return(list(looic = LOOIC, k_est = k_est))
    }

    }else{
      N <- dim(y)[1]

      if (ic=="dic"){

        k_vec <- c(1:k_max)
        DIC <- c()

        mod.mist.biv<-"model{
    # Likelihood:

    for (i in 1:N){
      yprev[i,1:2]<-y[i,1:2]
      y[i,1:2] ~ dmnorm(muOfClust[clust[i],],tauOfClust)
      clust[i] ~ dcat(pClust[1:k] )
    }

    # Prior:

    for (g in 1:k) {
      muOfClust[g,1:2] ~ dmnorm(mu_0[],S2[,])}
      tauOfClust[1:2,1:2] ~ dwish(S3[,],3)
      Sigma[1:2,1:2] <- inverse(tauOfClust[,])
      pClust[1:k] ~ ddirch(alpha)
  }"



        for (h in 2:k_max){
          clust_inits <- KMeans(y, k_vec[h])$cluster
          #cutree(hclust(dist(y), "average"),k)
          mu_inits <- matrix(0,k_vec[h],2)
          for (j in 1:k_vec[h]){
            mu_inits[j,] <- cbind(mean(y[clust_inits==j,1]), mean(y[clust_inits==j,2]))
          }
          #Reorder mu_inits according to the x-coordinate
          mu_inits <-
            mu_inits[sort(mu_inits[,1], decreasing=FALSE, index.return=TRUE)$ix,]
        # Data

        if (missing(priors)){
          mu_0 <- as.vector(c(0,0))
          S2 <- matrix(c(1,0,0,1),nrow=2)/100000
          S3 <- matrix(c(1,0,0,1),nrow=2)/100000
          alpha <- rep(1,k_vec[h])
        }else{
          if (is.null(priors$mu_0)){
            mu_0 <- as.vector(c(0,0))
          }else{
            mu_0 <- priors$mu_0
          }
          if (is.null(priors$S2)){
            S2 <- matrix(c(1,0,0,1),nrow=2)/100000
          }else{
            S2 <- priors$S2
          }
          if (is.null(priors$S3)){
            S3 <- matrix(c(1,0,0,1),nrow=2)/100000
          }else{
            S3 <- priors$S3
          }
          if (is.null(priors$alpha)){
            alpha <- rep(1, k_vec[h])
          }else{
            alpha <- priors$alpha
          }
        }

        init1.biv <- dump.format(list(muOfClust=mu_inits,
                                      tauOfClust= matrix(c(15,0,0,15),ncol=2),
                                      pClust=rep(1/k_vec[h],k_vec[h]), clust=clust_inits))
        init2.biv <- dump.format(list(muOfClust=mu_inits,
                                      tauOfClust= matrix(c(15,0,0,15),ncol=2),
                                      pClust=rep(1/k_vec[h],k_vec[h]), clust=clust_inits))
        init3.biv <- dump.format(list(muOfClust=mu_inits,
                                      tauOfClust= matrix(c(15,0,0,15),ncol=2),
                                      pClust=rep(1/k_vec[h],k_vec[h]), clust=clust_inits))
        dati.biv <- list(y = y, N = N, k = k_vec[h],
                         S2= S2, S3= S3, mu_0=mu_0,
                         alpha = alpha)
        moni.biv <- c("clust","muOfClust","tauOfClust","pClust", "dic")

        mod   <- mod.mist.biv
        dati  <- dati.biv
        init1 <- list(init1.biv, init2.biv, init3.biv)
        moni  <- moni.biv

        # Jags execution
        ogg.jags <- run.jags(model=mod, data=dati, monitor=moni,
                             inits=init1, n.chains=chains,plots=FALSE, thin=1,
                             sample=nMC, burnin=burn)
        DIC[h-1] <- as.double(ogg.jags$dic$dic)
        }
        k_est <- which.min(DIC)+1
        return(list(dic = DIC, k_est = k_est))


      }else if (ic == "loo"){
        if (missing(priors)){
          mu_0 <- c(0,0)
          eta <- 1
          sigma_d <- 2.5
        }else{
          if (is.null(priors$mu_0)){
            mu_0 <- c(0,0)
          }else{
            mu_0 <- priors$mu_0
          }
          if (is.null(priors$eta)){
            eta <- 1
          }else{
            eta  <- priors$eta
          }
          if (is.null(priors$sigma_d)){
            sigma_d <- 2.5
          }else{
            sigma_d <- priors$sigma_d
          }
        }

        k_vec <- c(1:k_max)
        LOOIC <- c()

        mix_biv <- "
        data {
          int<lower=1> k;          // number of mixture components
          int<lower=1> N;          // number of data points
          int D;                   // data dimension
          matrix[N,D] y;           // observations matrix
          vector[D] mu_0;
          real<lower=0> eta;
          real<lower=0> sigma_d;
        }
        parameters {
          simplex[k] theta;        // mixing proportions
          vector[D] mu[k];        // locations of mixture components
          cholesky_factor_corr[D] L_Omega;   // scales of mixture components
          vector<lower=0>[D] L_sigma;
          cholesky_factor_corr[D] L_tau_Omega;   // scales of mixture components
          vector<lower=0>[D] L_tau;
          }
        transformed parameters{
          vector[k] log_theta = log(theta);  // cache log calculation
          vector[k] pz[N];
          simplex[k] exp_pz[N];
          matrix[D,D] L_Sigma=diag_pre_multiply(L_sigma, L_Omega);
          matrix[D,D] L_Tau=diag_pre_multiply(L_tau, L_tau_Omega);


            for (n in 1:N){
                pz[n]=   multi_normal_cholesky_lpdf(y[n]|mu, L_Sigma)+
                         log_theta-
                         log_sum_exp(multi_normal_cholesky_lpdf(y[n]|
                                                     mu, L_Sigma)+
                         log_theta);
                exp_pz[n] = exp(pz[n]);
              }
          }
        model{
          L_Omega ~ lkj_corr_cholesky(eta);
          L_sigma ~ cauchy(0, sigma_d);
          mu ~ multi_normal_cholesky(mu_0, L_Tau);
            for (n in 1:N) {
              vector[k] lps = log_theta;
                for (j in 1:k){
                    lps[j] += multi_normal_cholesky_lpdf(y[n] |
                                                   mu[j], L_Sigma);
                    target+=pz[n,j];
                }
              target += log_sum_exp(lps);
          }
          }
        generated quantities{
          matrix[N,k] ll;
          vector[N] log_lik;
          int<lower=1, upper=k> z[N];

              for (n in 1:N){
                  z[n] = categorical_rng(exp_pz[n]);
                  for (j in 1:k){
                      ll[n,j] = multi_normal_cholesky_lpdf(y[n] | mu[j],
                                              L_Sigma)+log_theta[j];
              }
                log_lik[n] = log_sum_exp(ll[n,]);
                }
          }
          "


        for (h in 2:k_max){
          data =list(N=N, k=k_vec[h], y=y, D=2, mu_0=mu_0,
                     eta = eta, sigma_d = sigma_d)
          fit_biv <-  stan(model_code = mix_biv,
                           data=data,
                           chains =chains,
                           iter =nMC)
          sims_biv <- rstan::extract(fit_biv)
          log_lik <- extract_log_lik(fit_biv)
          loo <- loo(log_lik)
          LOOIC[h-1] <- loo$estimates[3,1]

        }
        k_est <- which.min(LOOIC)+1
        return(list(looic=LOOIC, k_est = k_est))
      }

  }


}

