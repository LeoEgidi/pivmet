## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load, warning =FALSE, message = FALSE-------------------------------
library(pivmet)

## ----fish_hist, fig.align ='center'--------------------------------------
data(fish)
y <- fish[,1]
hist(y, breaks=40, prob = TRUE, cex.lab=1.6,
            main ="Fishery data", cex.main =1.7,
            col="navajowhite1", border="navajowhite1")
lines(density(y), lty=1, lwd=3, col="blue")

## ----fish_data-----------------------------------------------------------
library(pivmet)
data(fish)
y <- fish[,1]
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

## ----nested, fig.align ='center'-----------------------------------------
set.seed(50)
N  <- 200
k  <- 3
nMC <- 2000
M1 <- c(-10,8)
M2 <- c(10,.1)
M3 <- c(30,8)
Mu <- matrix(rbind(M1,M2,M3),c(k,2))
stdev    <- cbind(rep(1,k), rep(200,k))
Sigma.p1 <- matrix(c(stdev[1,1],0,0,stdev[1,1]),
nrow=2, ncol=2)
Sigma.p2 <- matrix(c(stdev[1,2],0,0,stdev[1,2]),
 nrow=2, ncol=2)
W   <- c(0.2,0.8)
sim <- piv_sim(N = N, k = k, Mu = Mu,
  stdev = stdev, Sigma.p1 = Sigma.p1, Sigma.p2 = Sigma.p2, W = W)
plot(sim$y, xlab="y[,1]", ylab="y[,2]")

## ----mcmc, message =  FALSE, warning = FALSE-----------------------------
res <- piv_MCMC(y = sim$y, k= k, nMC =nMC)

## ----pivotal_rel, fig.show='hold', fig.align='center',fig.width=7--------
rel <- piv_rel(mcmc=res)
piv_plot(y = sim$y, mcmc = res, rel_est = rel, par = "mean", type = "chains")
piv_plot(y = sim$y, mcmc = res, rel_est = rel, type = "hist")

## ----fish_data-----------------------------------------------------------
library(pivmet)
data(fish)
y <- fish[,1]
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

## ----nested, fig.align ='center'-----------------------------------------
set.seed(50)
N  <- 200
k  <- 3
nMC <- 2000
M1 <- c(-10,8)
M2 <- c(10,.1)
M3 <- c(30,8)
Mu <- matrix(rbind(M1,M2,M3),c(k,2))
stdev    <- cbind(rep(1,k), rep(200,k))
Sigma.p1 <- matrix(c(stdev[1,1],0,0,stdev[1,1]),
nrow=2, ncol=2)
Sigma.p2 <- matrix(c(stdev[1,2],0,0,stdev[1,2]),
 nrow=2, ncol=2)
W   <- c(0.2,0.8)
sim <- piv_sim(N = N, k = k, Mu = Mu,
  stdev = stdev, Sigma.p1 = Sigma.p1, Sigma.p2 = Sigma.p2, W = W)
plot(sim$y, xlab="y[,1]", ylab="y[,2]")

## ----mcmc, message =  FALSE, warning = FALSE-----------------------------
res <- piv_MCMC(y = sim$y, k= k, nMC =nMC)

## ----pivotal_rel, fig.show='hold', fig.align='center',fig.width=7--------
rel <- piv_rel(mcmc=res)
piv_plot(y = sim$y, mcmc = res, rel_est = rel, par = "mean", type = "chains")
piv_plot(y = sim$y, mcmc = res, rel_est = rel, type = "hist")

