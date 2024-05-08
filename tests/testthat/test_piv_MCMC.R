context("arguments errors")
test_that("piv_MCMC recognise errors/warnings",{

   N   <- 200
   k   <- 4
   D   <- 2
   nMC <- 1000
   M1  <- c(-.5,8)
   M2  <- c(25.5,.1)
   M3  <- c(49.5,8)
   M4  <- c(63.0,.1)
   Mu  <- rbind(M1,M2,M3,M4)
   Sigma.p1 <- diag(D)
   Sigma.p2 <- 20*diag(D)
   W <- c(0.2,0.8)
   sim <- piv_sim(N = N, k = k, Mu = Mu,
                  Sigma.p1 = Sigma.p1,
                  Sigma.p2 = Sigma.p2, W = W)
   # k = 0
   expect_error(piv_MCMC(y = sim$y, k =0, nMC = nMC))
   # wrong software
   expect_error(piv_MCMC(y = sim$y, k =k, nMC = nMC, software = "some"))
   # wrong piv_criterion
   expect_error(piv_MCMC(y = sim$y, k =k, nMC = nMC, piv.criterion = "some"))
   # wrong clustering method
   expect_error(piv_MCMC(y = sim$y, k =k, nMC = nMC, clustering ="kmeans"))
   # wrong prior choice
   expect_error(piv_MCMC(y = sim$y, k =k, nMC = nMC, priors = "normal"))
   # wrong prior choice dimensions (for mu_0)
   expect_error(piv_MCMC(y = sim$y, k =k, nMC = nMC, priors = list(mu_0 = c(1,2,3) )))
   # wrong prior choice dimensions (for alpha)
   expect_error(piv_MCMC(y = sim$y, k =k, nMC = nMC, priors = list(mu_0 = c(1,2), alpha = rep(1, k+1) )))
   # wrong prior choice + software
   expect_error(piv_MCMC(y = sim$y, k =k, nMC = nMC,
                         priors = list(kind = "independence", parameter = "priorsFish", hierarchical = "tau"),
                          sofware = "rstan"))


})
