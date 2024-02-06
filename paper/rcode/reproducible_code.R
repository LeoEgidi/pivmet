# load pivmet package

library(pivmet)

############################################
# Example 1: relabeling on fishery data
############################################

library(bayesmix)
set.seed(100)

# load data
data(fish)
y <- fish[,1]
k <- 5
nMC <- 15000

#fit the mixture model for univariate data and find the pivots
res <- piv_MCMC(y = y, k = k, nMC = nMC, burn = 0.5*nMC, software = "rjags")

# relabel the chains: figure 2
rel <- piv_rel(mcmc=res)
piv_plot(y = y, mcmc = res, rel_est = rel, type="chains")

# use Stan
res_stan <- piv_MCMC(y = y, k = k, nMC = nMC/3, burn = 0.5*nMC/3, software ="rstan")
cat(res_stan$model)



#######################################################
# Example 2: consensus clustering
#######################################################

library(mclust)
library(cluster)
library(mvtnorm)


# simulate data
set.seed(123)
n=620
k=3
n1=20
n2=100
n3=500
x=matrix(NA, n,2)
gruppovero=c( rep(1,n1), rep(2, n2), rep(3, n3))

for (i in 1:n1){
  x[i,]=rmvnorm(1, c(1,5), sigma=diag(2))
}
for (i in 1:n2){
  x[n1+i,]=rmvnorm(1, c(4,0), sigma=diag(2))
}
for (i in 1:n3){
  x[n1+n2+i,]=rmvnorm(1, c(6,6), sigma=diag(2))
}

kmeans_res <- kmeans(x, centers=k)

res1 <- piv_KMeans(x, k, alg.type = "hclust",
                   piv.criterion ="maxsumdiff",
                   prec_par=n1)

pivots1 <- res1$pivots

res2 <- piv_KMeans(x, k, alg.type = "hclust",
                   piv.criterion ="maxsumint",
                   prec_par=n1)

pivots2 <- res2$pivots

res3 <- piv_KMeans(x, k, alg.type = "hclust",
                   piv.criterion ="minsumnoint",
                   prec_par = n1)

pivots3 <- res3$pivots

res4 <- piv_KMeans(x, k,
                   alg.type = "hclust",
                   piv.criterion ="MUS",
                   prec_par=n1+n2)

pivots4 <- res4$pivots

a.wgt  <- agnes(x, method = "average")
a.sing <- agnes(x, method = "single")
a.comp <- agnes(x, method = "complete")
pamx <- pam(x, k)

partition=list(kmeans_res$cluster, res1$cluster,  res2$cluster,
               res4$cluster,  pamx$clustering, cutree(a.wgt, k),
               cutree(a.sing, k), cutree(a.comp, k))

adj <- c()
for(l in 1: length(partition))
{
  adj[l]=adjustedRandIndex(gruppovero, partition[[l]])
}

a<-round(adj,3)



# Figure 3

par(mfrow=c(3,3), oma=c(0.1,0.1,0.1,0.1), mar=c(5,4.3,3,1), pty="s")
colors_cluster=c("grey", "darkolivegreen3", "coral")
colors_centers=c("black", "darkgreen", "firebrick")
plot(x, col = colors_cluster[gruppovero]
     ,bg= colors_cluster[gruppovero], pch=21, xlab="y[,1]",
     ylab="y[,2]", cex.lab=1.5,
     main="(a) True data", cex.main=1.5)
plot(x, col = colors_cluster[kmeans_res$cluster],
     bg=colors_cluster[kmeans_res$cluster], pch=21, xlab="y[,1]", ylab="y[,2]",
     cex.lab=1.5,main="(b) K-means",  cex.main=1.5)
points(kmeans_res$centers, col = colors_centers[1:k], pch = 8, cex = 2)

plot(x, col = colors_cluster[res4$cluster], bg=colors_cluster[res4$cluster], pch=21, xlab="y[,1]", ylab="y[,2]", cex.lab=1.5,
     main="(c) piv_KMeans[1]", cex.main=1.5)
points(x[pivots4[1],1], x[pivots4[1],2], pch=24, col=colors_centers[1],bg=colors_centers[1],   cex=1.5)
points(x[pivots4[2],1], x[pivots4[2],2], pch=24,  col=colors_centers[2], bg=colors_centers[2], cex=1.5)
points(x[pivots4[3],1], x[pivots4[3],2], pch=24, col=colors_centers[3], bg=colors_centers[3], cex=1.5)
points(res4$centers, col = colors_centers[1:k], pch = 8, cex = 2)


plot(x, col = colors_cluster[res2$cluster], bg=colors_cluster[res2$cluster], pch=21, xlab="y[,1]", ylab="y[,2]", cex.lab=1.5,
     main="(d) piv_KMeans[2]", cex.main=1.5)
points(x[pivots2[1],1], x[pivots2[1],2], pch=24, col=colors_centers[1],bg=colors_centers[1],   cex=1.5)
points(x[pivots2[2],1], x[pivots2[2],2], pch=24,  col=colors_centers[2], bg=colors_centers[2], cex=1.5)
points(x[pivots2[3],1], x[pivots2[3],2], pch=24, col=colors_centers[3], bg=colors_centers[3], cex=1.5)
points(res2$centers, col = colors_centers[1:k], pch = 8, cex = 2)


plot(x, col = colors_cluster[res1$cluster], bg=colors_cluster[res1$cluster], pch=21, xlab="y[,1]", ylab="y[,2]", cex.lab=1.5,
     main="(e) piv_KMeans[4]", cex.main=1.5)
points(x[pivots1[1],1], x[pivots1[1],2], pch=24, col=colors_centers[1],bg=colors_centers[1],   cex=1.5)
points(x[pivots1[2],1], x[pivots1[2],2], pch=24,  col=colors_centers[2], bg=colors_centers[2], cex=1.5)
points(x[pivots1[3],1], x[pivots1[3],2], pch=24, col=colors_centers[3], bg=colors_centers[3], cex=1.5)
points(res1$centers, col = colors_centers[1:k], pch = 8, cex = 2)


plot(x, col=colors_cluster[pamx$clustering], bg =colors_cluster[pamx$clustering],
     main="(g) PAM",
     pch=21, xlab="y[,1]", ylab="y[,2]", cex.lab=1.5,
     cex.main=1.5)

plot(x, col=colors_cluster[cutree(a.wgt, 3)], bg = colors_cluster[cutree(a.wgt, 3)],
     main="(h) Agnes - average", pch=21, xlab="y[,1]", ylab="y[,2]", cex.lab=1.5,
     cex.main=1.5)


plot(x, col=colors_cluster[cutree(a.sing, 3)], bg = colors_cluster[cutree(a.sing, 3)],
     main="(i) Agnes - single", pch=21, xlab="y[,1]", ylab="y[,2]", cex.lab=1.5,
     cex.main=1.5)

plot(x, col=colors_cluster[cutree(a.comp, 3)], bg = colors_cluster[cutree(a.comp, 3)],
     main="(l) Agnes - complete", pch=21, xlab="y[,1]", ylab="y[,2]", cex.lab=1.5,
     cex.main=1.5)








