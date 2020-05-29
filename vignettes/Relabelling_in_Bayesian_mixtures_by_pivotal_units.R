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

