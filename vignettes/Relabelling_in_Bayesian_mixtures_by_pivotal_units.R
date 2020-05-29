## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load, warning =FALSE, message = FALSE-------------------------------
library(pivmet)
library(bayesmix)

## ----nested, fig.align ='center'-----------------------------------------
set.seed(50)
N  <- 200
k  <- 3
nMC <- 2000
M1 <- c(-10,8)
M2 <- c(10,.1)
M3 <- c(30,8)
# matrix of input means
Mu <- matrix(rbind(M1,M2,M3),c(k,2)) 
sds    <- cbind(rep(1,k), rep(15,k))
# covariance matrices for the two subgroups
Sigma.p1 <- matrix(c(sds[1,1]^2,0,0,sds[1,1]^2),
nrow=2, ncol=2)
Sigma.p2 <- matrix(c(sds[1,2]^2,0,0,sds[1,2]^2),
nrow=2, ncol=2)
# subgroups' weights
W   <- c(0.2,0.8)
# simulate data
sim <- piv_sim(N = N, k = k, Mu = Mu,
 Sigma.p1 = Sigma.p1, Sigma.p2 = Sigma.p2, W = W)

## ----mcmc, message =  FALSE, warning = FALSE-----------------------------
res <- piv_MCMC(y = sim$y, k= k, nMC =nMC, 
                piv.criterion = "maxsumdiff")

## ----pivotal_rel, fig.show='hold', fig.align='center',fig.width=7--------
rel <- piv_rel(mcmc=res)
piv_plot(y = sim$y, mcmc = res, rel_est = rel, par = "mean", type = "chains")
piv_plot(y = sim$y, mcmc = res, rel_est = rel, type = "hist")

## ----fish_hist, fig.align ='center', fig.width=5.5-----------------------
data(fish)
y <- fish[,1]
hist(y, breaks=40, prob = TRUE, cex.lab=1.6,
             main ="Fishery data", cex.main =1.7,
             col="navajowhite1", border="navajowhite1")
 lines(density(y), lty=1, lwd=3, col="blue")

## ----fish_data-----------------------------------------------------------
k <- 5
nMC <- 15000
res <- piv_MCMC(y = y, k = k, nMC = nMC, burn = 0.5*nMC,
software = "rjags")

## ----true_iter-----------------------------------------------------------
res$true.iter

## ----fish_rel, fig.align= 'center', fig.width=7--------------------------
rel <- piv_rel(mcmc=res)
piv_plot(y=y, res, rel, par = "mean", type="chains")
piv_plot(y=y, res, rel, type="hist")

## ----model_code----------------------------------------------------------
cat(res$model)

