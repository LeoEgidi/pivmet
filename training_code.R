
## Training code to be launched for a brief package illustration.
## The following code is protected by the GPL-2.0 License



## packages installation and loading-------------------------------------------------------------------
install.packages("devtools")
install.packages("bayesmix")
install.packages("mvtnorm")
library(bayesmix)
library(devtools)
library(mvtnorm)
devtools::install_github("leoegidi/pivmet")  # github installation

## alternative: CRAN installation, containing an older version of the package
# install.packages("pivmet)


## ----example 1: mixture modeling with fish dataset-----------------------------------------------------------------------------------------
library(pivmet)
data(fish)
y <- fish[,1]
N <- length(y)  # sample size
k <- 5          # fixed number of clusters
nMC <- 12000    # MCMC iterations


## ----MCMC fit using the underlying JAGS sampling------------------------------------------------------------
res <- piv_MCMC(y = y, k = k, nMC = nMC)


## ----relabeling step-----------------------------------------------------------
rel <- piv_rel(mcmc=res)
piv_plot(y = y, mcmc = res, rel_est = rel, type = "chains")
piv_plot(y = y, mcmc = res, rel_est = rel, type = "hist")


## ----sparse finite mixtures-------------------------------------------------------
res2 <- piv_MCMC(y, k, nMC, sparsity = TRUE,
                 priors = list(alpha = rep(0.001, k))) # sparse on eta
barplot(table(res2$nclusters), xlab= expression(K["+"]),
        col = "blue", border = "red", main = expression(paste("p(",K["+"], "|y)")),
        cex.main=3, yaxt ="n", cex.axis=2.4, cex.names=2.4,
        cex.lab=2)


## ----example 2: k-means with simulated data----------------------------------

set.seed(123)
n  <- 620
centers  <- 3
n1 <- 20
n2 <- 100
n3 <- 500
x  <- matrix(NA, n,2)
truegroup <- c( rep(1,n1), rep(2, n2), rep(3, n3))

for (i in 1:n1){
 x[i,]=rmvnorm(1, c(1,5), sigma=diag(2))}
for (i in 1:n2){
 x[n1+i,]=rmvnorm(1, c(4,0), sigma=diag(2))}
for (i in 1:n3){
 x[n1+n2+i,]=rmvnorm(1, c(6,6), sigma=diag(2))}

H <- 1000
a <- matrix(NA, H, n)

  for (h in 1:H){
    a[h,] <- kmeans(x,centers)$cluster
  }

#build the similarity matrix
sim_matr <- matrix(NA, n,n)
 for (i in 1:(n-1)){
    for (j in (i+1):n){
      sim_matr[i,j] <- sum(a[,i]==a[,j])/H
      sim_matr[j,i] <- sim_matr[i,j]
    }
  }

cl <- kmeans(x, centers, nstart=10)$cluster
mus_alg <- MUS(C = sim_matr, clusters = cl, prec_par = 5)


## ----kmeans_plots--------------------------------------------
# launch classical kmeans
kmeans_res <- kmeans(x, centers, nstart = 10)
# plots
par(mfrow=c(1,2))
colors_cluster <- c("grey", "darkolivegreen3", "coral")
colors_centers <- c("black", "darkgreen", "firebrick")

graphics::plot(x, col = colors_cluster[truegroup]
                 ,bg= colors_cluster[truegroup], pch=21,
                  xlab="y[,1]",
                  ylab="y[,2]", cex.lab=1.5,
                  main="True data", cex.main=1.5)

graphics::plot(x, col = colors_cluster[kmeans_res$cluster],
      bg=colors_cluster[kmeans_res$cluster], pch=21, xlab="y[,1]",
      ylab="y[,2]", cex.lab=1.5,main="K-means",  cex.main=1.5)
points(kmeans_res$centers, col = colors_centers[1:centers],
      pch = 8, cex = 2)


## ----piv_KMeans: robust k-means with pivotal seeding---------------------------------------------------------------------------
piv_res <- piv_KMeans(x, centers)
# plots
par(mfrow=c(1,2), pty="s")
colors_cluster <- c("grey", "darkolivegreen3", "coral")
colors_centers <- c("black", "darkgreen", "firebrick")
graphics::plot(x, col = colors_cluster[truegroup],
   bg= colors_cluster[truegroup], pch=21, xlab="x[,1]",
   ylab="x[,2]", cex.lab=1.5,
   main="True data", cex.main=1.5)

graphics::plot(x, col = colors_cluster[piv_res$cluster],
   bg=colors_cluster[piv_res$cluster], pch=21, xlab="x[,1]",
   ylab="x[,2]", cex.lab=1.5,
   main="piv_Kmeans", cex.main=1.5)
points(x[piv_res$pivots[1],1], x[piv_res$pivots[1],2],
   pch=24, col=colors_centers[1],bg=colors_centers[1],
   cex=1.5)
points(x[piv_res$pivots[2],1], x[piv_res$pivots[2],2],
   pch=24,  col=colors_centers[2], bg=colors_centers[2],
   cex=1.5)
points(x[piv_res$pivots[3],1], x[piv_res$pivots[3],2],
   pch=24, col=colors_centers[3], bg=colors_centers[3],
   cex=1.5)
points(piv_res$centers, col = colors_centers[1:centers],
   pch = 8, cex = 2)


