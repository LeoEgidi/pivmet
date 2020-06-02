## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load, warning =FALSE, message = FALSE------------------------------------
library(pivmet)
library(bayesmix)
library(bayesplot)

## ----nested, fig.align ='center'----------------------------------------------
set.seed(500)
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
Sigma.p2 <- 15*diag(D)
# subgroups' weights
W   <- c(0.2,0.8)
# simulate data
sim <- piv_sim(N = N, k = k, Mu = Mu,
 Sigma.p1 = Sigma.p1, Sigma.p2 = Sigma.p2, W = W)

