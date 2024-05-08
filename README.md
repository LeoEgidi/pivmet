
<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
)
```

# pivmet

The goal of ```pivmet``` is to propose some pivotal methods in order to:

- undo the label switching problem which naturally arises during the MCMC sampling in Bayesian mixture models $\rightarrow$ **pivotal relabelling** (Egidi et al. 2018a)

- fit sparse finite Gaussian mixtures

- initialize the K-means algorithm aimed at obtaining a good clustering solution $\rightarrow$ **pivotal seeding** (Egidi et al. 2018b)

## Installation

- <span style="color:red">PAY ATTENTION! BEFORE INSTALLING</span>: make sure to download the JAGS program at
[https://sourceforge.net/projects/mcmc-jags/](https://sourceforge.net/projects/mcmc-jags/).

You can install the CRAN version of ```pivmet``` with:

```{r, eval = FALSE}
install.packages("pivmet")
library(pivmet)
```

You can install the development version of ```pivmet```  from Github with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("leoegidi/pivmet")
```

## Example 1. Dealing with label switching: relabelling in Bayesian mixture models by pivotal units (fish data)

First of all, we load the package and we import the ```fish``` dataset belonging to the ```bayesmix``` package:

```{r example}
library(bayesmix)
library(pivmet)
data(fish)
y <- fish[,1]
N <- length(y)  # sample size 
k <- 5          # fixed number of clusters
nMC <- 12000    # MCMC iterations
```

Then we fit a Bayesian Gaussian mixture using the ```piv_MCMC``` function:


```{r fit, message =FALSE, warning = FALSE}
res <- piv_MCMC(y = y, k = k, nMC = nMC)
```


Finally, we can apply pivotal relabelling and inspect the new posterior estimates with the functions ```piv_rel``` and ```piv_plot```, respectively:

```{r plot, message =FALSE, warning = FALSE}
rel <- piv_rel(mcmc=res)
piv_plot(y = y, mcmc = res, rel_est = rel, type = "chains")
piv_plot(y = y, mcmc = res, rel_est = rel, type = "hist")
```


To allow sparse finite mixture fit, we could select the argument ```sparsity = TRUE```:


```{r sparsity, message =FALSE, warning = FALSE}
res2 <- piv_MCMC(y, k, nMC, sparsity = TRUE,
                 priors = list(alpha = rep(0.001, k))) # sparse on eta
barplot(table(res2$nclusters), xlab= expression(K["+"]),
        col = "blue", border = "red", main = expression(paste("p(",K["+"], "|y)")),
        cex.main=3, yaxt ="n", cex.axis=2.4, cex.names=2.4,
        cex.lab=2)
```



## Example 2. K-means clustering using MUS and other pivotal algorithms

Sometimes K-means algorithm does not provide an optimal clustering solution. Suppose to generate some clustered data and to detect one pivotal unit for each group with the ```MUS``` (Maxima Units Search algorithm) function:

```{r mus, echo =TRUE, eval = TRUE, message = FALSE, warning = FALSE}
library(mvtnorm)

#generate some data

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
```


Quite often, classical K-means fails in recognizing the *true* groups:



```{r kmeans_plots, echo =TRUE, fig.show='hold', eval = TRUE, message = FALSE, warning = FALSE}
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
```


In such situations, we may need a more robust version of the classical K-means. The pivots may be used as initial seeds for a classical K-means algorithm. The function `piv_KMeans` works as the classical `kmeans` function, with some optional arguments (in the figure below, the colored triangles represent the pivots).

```{r musk, fig.show='hold'}
# launch piv_KMeans
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

```

## References

Egidi, L., Pappadà, R., Pauli, F. and Torelli, N. (2018a). Relabelling in Bayesian Mixture Models by Pivotal Units. Statistics and Computing, 28(4), 957-969.

Egidi, L., Pappadà, R., Pauli, F., Torelli, N. (2018b). K-means seeding via MUS algorithm. Conference Paper, Book of Short Papers, SIS2018, ISBN: 9788891910233.
