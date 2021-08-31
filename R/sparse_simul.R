#
N   <- 100
k   <- 4
k_over <- 15
D   <- 4
nMC <- 2000
M1  <- c(2,-2,0,0)
M2  <- -M1
M3  <- c(2,2,0,0)
M4  <- -M3
Mu  <- rbind(M1,M2,M3,M4)
Sigma.p1 <- diag(D)
Sigma.p2 <- diag(D)
W <- c(0.5,0.5)
sim <- piv_sim(N = N, k = k, Mu = Mu,
               Sigma.p1 = Sigma.p1,
               Sigma.p2 = Sigma.p2, W = W)
res <- piv_MCMC(y = sim$y, k = k_over, nMC = nMC,
                sparsity =TRUE,
                priors = list(a_sp=2, b_sp=2))

# darwin data (bayesmix)
library(bayesmix)
data(darwin)
y <- darwin$darwin
nMC <- 5000
k <- 8


## rjags
res2 <- piv_MCMC(y, k, nMC, sparsity = TRUE,
                 priors = list(alpha = rep(0.01, k))) # sparse on eta
pdf(file="sparse_001.pdf", height = 6, width = 8)
barplot(table(res2$nclusters), xlab= expression(K["+"]),
        col = "blue", border = "red", main = expression(paste("p(",K["+"], "|y)")),
        cex.main=2, yaxt ="n", cex.axis=1.6, cex.names=1.6,
        cex.lab=2)
dev.off()


res <- piv_MCMC(y, k, nMC, sparsity = TRUE) # uniform on eta
pdf(file="sparse_1.pdf", height = 6, width = 8)
barplot(table(res$nclusters), xlab= expression(K["+"]),
        col = "blue", border = "red", main = expression(paste("p(",K["+"], "|y)")),
        cex.main=2, yaxt ="n",cex.axis=1.6, cex.names=1.6,
        cex.lab=2)
dev.off()

res3 <- piv_MCMC(y, k, nMC, sparsity = TRUE,
                 priors = list(alpha = rep(10, k))) # overfitting on eta
pdf(file="sparse_10.pdf", height = 6, width = 8)
barplot(table(res3$nclusters), xlab= expression(K["+"]),
        col = "blue", border = "red", main = expression(paste("p(",K["+"], "|y)")),
        cex.main=2, yaxt ="n", cex.axis=1.6, cex.names=1.6,
        cex.lab=2)
dev.off()


