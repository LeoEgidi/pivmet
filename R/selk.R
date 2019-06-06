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

piv_k <- function(y, priors,  software =c("rjags", "rstan"),
                  nMC = 500, chains = 4, kmax = 8){
  if (is.vector(y)){
    if (software=="rjags"){
      if (missing(priors)){
        priors = list(kind = "independence",
                      parameter = "priorsFish",
                      hierarchical = "tau")
      }

    k_vec <- c(1:k_max)
    DIC <- c()
    for (h in 2:kmax){
    mod.mist.univ <- BMMmodel(y, k = k_vec[h],
                              initialValues = list(S0 = 2),
                              priors = priors)
    control <- JAGScontrol(variables = c("mu", "tau", "eta", "S", "dic"),
                           burn.in = burn, n.iter = nMC, seed = 10)
    ogg.jags <- JAGSrun(y, model = mod.mist.univ, control = control)

    }
    }else if (software == "rstan"){
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

      for (h in 5:k_max){
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

    }
    }else{

  }

       return(looic = LOOIC)
}
