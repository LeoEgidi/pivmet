## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load, warning =FALSE, message = FALSE-------------------------------
library(pivmet)
library(bayesmix)
library(bayesplot)

## ----nested, fig.align ='center'-----------------------------------------
set.seed(50)
N  <- 200
k  <- 3
D <- 2
nMC <- 2000
M1 <- c(-10,8)
M2 <- c(10,.1)
M3 <- c(30,8)
# matrix of input means
Mu <- rbind(M1,M2,M3)
# covariance matrices for the two subgroups
Sigma.p1 <- diag(D)
Sigma.p2 <- (15^2)*diag(D)
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
res <- piv_MCMC(y = y, k = k, nMC = nMC, 
                burn = 0.5*nMC, software = "rjags")

## ----true_iter-----------------------------------------------------------
res$true.iter

## ----fish_rel, fig.align= 'center', fig.width=7--------------------------
rel <- piv_rel(mcmc=res)
piv_plot(y=y, res, rel, par = "mean", type="chains")
piv_plot(y=y, res, rel, type="hist")

## ----stan, fig.align= 'center', fig.width=7------------------------------
res2 <- piv_MCMC(y = y, k = k, nMC = 3000, 
                 software = "rstan")
rel2 <- piv_rel(res2)
piv_plot(y=y, res2, rel2, par = "mean", type="chains")

## ----bayesplot,fig.align= 'center', fig.width=5--------------------------
posterior <- as.array(res2$stanfit)
mcmc_intervals(posterior, regex_pars = c("mu"))

## ----model_code----------------------------------------------------------
cat(res$model)
cat(res2$model)

