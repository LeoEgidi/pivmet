# load pivmet package

library(pivmet)

#####################################################
## 3.3. Mixtures of bivariate Gaussian distributions
#####################################################

library(rstan)
set.seed(10)
# Simulate data
N <- 150
k <- 4
D <- 2
nMC <- 5000
M1  <- c(-.5,8)
M2  <- c(25.5,.1)
M3  <- c(49.5,8)
M4  <- c(25,25)
Mu  <- rbind(M1,M2,M3,M4)
Sigma.p1 <- diag(D)
Sigma.p2 <- (14^2)*diag(D)
W <- c(0.2,0.8)

#Simulate data
sim <- piv_sim(N = N, k = k , Mu = Mu, Sigma.p1 = Sigma.p1, Sigma.p2 = Sigma.p2, W = W)

# Pivotal relabelling
res <- piv_MCMC(y = sim$y, k = k, nMC =nMC, piv.criterion="MUS", software="rstan")
rel <- piv_rel(mcmc = res)

# Plot the Results
piv_plot(y=sim$y, mcmc=res, rel_est = rel, par ="mean", type="chains")
piv_plot(y=sim$y, mcmc=res, rel_est = rel, type="hist")


#######################################################
#3.4. Consensus clustering based on pivots
#######################################################

library(foreign)
library(mclust)

# read data
file.n= "https://raw.githubusercontent.com/deric/clustering-benchmark/master/src/main/resources/datasets/artificial/2d-3c-no123.arff"
data=read.arff(file.n)
x=data[, 1:2]

# Perform pivKmeans on x with 3 groups
pkm=piv_KMeans(x, 3, alg.type = "hclust")
layout(1:1)
# visualize the clusters
plot(x, col = pkm$cluster, pch=19)
legend("bottomright", legend=c("1", "2", "3"), pch=19, col=c(1:3))
#Compare with 'true' labels
table(pkm$cluster, as.numeric(data[, 3]))
adjustedRandIndex(pkm$cluster, as.numeric(data[, 3]))

#########################################
# 3.5. Dirichlet mixture process
########################################

library(dirichletprocess)
library(ggplot2)
set.seed(1234)
n <- 200 # sample size
nMC <- 1000 # MCMC iterations
y <- rt(n, 3) + 2 #generate sample data
dp <- DirichletProcessGaussian(y) 
dp <- Fit(dp, nMC) # MCMC sampling (dirichletprocess)
dp$numberClusters # number of "non-empty" clusters

C_array <- array(1, dim = c(nMC, n, n)) # co-association array
for (h in 1:nMC){
  for (i in 1:(n-1)){
    for (j in (i+1):n){
      if (dp$labelsChain[[h]][i]==dp$labelsChain[[h]][j]){
      C_array[h,i,j] <- 1
      }else{
      C_array[h,i,j] <- 0  
      }
    }
  }
}

C <- apply(C_array, c(2,3), mean) # co-association matrix 
piv_selection <- piv_sel(C = C, clusters = dp$clusterLabels) # pivotal methods
piv_index <- piv_selection$pivots[,3] # maxsumdiff
df_piv <- data.frame(x=y[piv_index], y = rep(-0.01, dp$numberClusters))
plot(dp)+geom_point(aes(x = x, y = y), data = df_piv,
                    colour = "blue", size = 2.5)
ggsave(file="dpmm_pivots.pdf", width =8, height =8)

##############################
# 4.1. Fishery Data
##############################
library(bayesmix)
set.seed(100)
# load data
data(fish)
y <- fish[,1]
k <- 5
nMC <- 15000

#fit the mixture model for univariate data and find the pivots
res <- piv_MCMC(y = y, k = k, nMC = nMC, burn = 0.5*nMC, software = "rjags")

# relabel the chains
rel <- piv_rel(mcmc=res)
piv_plot(y = y, mcmc = res, rel_est = rel, type="chains")

#piv_plot(y = y, mcmc = res, rel_est = rel, type="hist")

# use Stan
library(rstan)
res_stan <- piv_MCMC(y = y, k = k, nMC = nMC/3, burn = 0.5*nMC/3, software ="rstan")
cat(res_stan$model)


####################
# 4.2. Galaxies Data
####################

library(MASS)
set.seed(101)

#load and trasform data
data(galaxies)
y <- galaxies
y <- y/1000

#Pivotal relabelling

k <- 3
nMC <- 15000
res <- piv_MCMC(y = y, k = k, nMC = nMC, priors=list(B0inv = 10), burn = 0.5*nMC)
rel <- piv_rel(mcmc=res)
piv_plot(y=y, res, rel, type="chains")
