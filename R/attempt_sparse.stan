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
