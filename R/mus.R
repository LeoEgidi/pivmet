#'MUS algorithm
#'
#' Perform Maxima Units Search (MUS) algorithm on a large and sparse matrix in
#' order to find a set of pivotal units through a sequential search
#' in the given matrix.
#'
#' @param C  \eqn{N \times N} matrix with a non-negligible number of zeros.
#' For instance, a similarity matrix estimated from a \eqn{N \times D} data matrix whose rows
#' are statistical units, or a co-association matrix resulting from clustering
#' ensembles.
#' @param clusters A vector of integers from \eqn{1:k}
#' indicating the cluster to which each point is allocated (it requires \eqn{k < 5}, see Details).
#' @param prec_par  Optional argument. The maximum number of candidate pivots for each group.
#' Default is 10.
#' @details
#'
#' Consider \eqn{H} distinct partitions of a set of \eqn{N} \eqn{d}-dimensional
#' statistical units into \eqn{k}
#' groups determined by some
#' clustering technique.  A \eqn{N \times N} co-association matrix
#' \eqn{C} with generic element \eqn{c_{i,j}=n_{i,j}/H} can be constructed,
#' where \eqn{n_{i,j}} is the number of times the \eqn{i}-th and the \eqn{j}-th unit
#' are assigned to the same cluster with respect to the clustering ensemble.
#' Units which are very distant
#' from each other are likely to have zero co-occurrences; as a consequence,
#'  \eqn{C} is
#' a square symmetric matrix expected  to contain a non-negligible number of zeros.
#' The main task of the MUS algorithm is to detect submatrices of small
#' rank from the co-association matrix
#' and extract those units---pivots---such
#' that the \eqn{k \times k} submatrix of \eqn{C},
#' determined by only the pivotal rows
#' and columns indexes, is identical or nearly identical.
#' Practically, the resulting units
#' have the desirable property to be representative of
#' the group they belong to.
#'
#' With the argument \code{prec_par} the user may increase
#' the powerful of the underlying MUS algorithm (see @egidi2018mus for details).
#' Given the default value 10, the function internally computes an
#' effective \code{prec_par} as \eqn{\min( 10, \min n_j )},
#' where \eqn{n_j} is the number of units belonging to the group
#' \eqn{j, \ j=1,\ldots,k}.
#'
#' @return
#'
#' \item{\code{pivots}}{ A vector of integers in 1:N denoting the indeces of the \code{k} selcted pivotal units.}
#' \item{\code{prec_par}}{The effective number of alternative pivots considered for each group. See Details.}
#'
#' @author Leonardo Egidi \email{legidi@units.it}, Roberta Pappadà
#' @references Egidi, L., Pappadà, R., Pauli, F., Torelli, N. (2018).
#'  Maxima Units Search(MUS) algorithm:
#' methodology and applications. In: Perna, C. , Pratesi, M., Ruiz-Gazen A. (eds.) Studies in
#' Theoretical and Applied Statistics,
#' Springer Proceedings in Mathematics and Statistics 227, pp. 71–81.
#'
#'
#'
#' @examples
#'
#' # Data generated from a mixture of three bivariate Gaussian distributions
#'
#'\dontrun{
#' N <- 620
#' centers  <- 3
#' n1 <- 20
#' n2 <- 100
#' n3 <- 500
#' x  <- matrix(NA, N,2)
#' truegroup <- c( rep(1,n1), rep(2, n2), rep(3, n3))
#'
#'
#'  x[1:n1,]=rmvnorm(n1, c(1,5), sigma=diag(2))
#'  x[(n1+1):(n1+n2),]=rmvnorm(n2, c(4,0), sigma=diag(2))
#'  x[(n1+n2+1):(n1+n2+n3),]=rmvnorm(n3, c(6,6), sigma=diag(2))
#'
#' # Build a similarity matrix from clustering ensembles
#'
#' H <- 1000
#' a <- matrix(NA, H, N)
#'
#' for (h in 1:H){
#'    a[h,] <- kmeans(x,centers)$cluster
#' }
#'
#' sim_matr <- matrix(NA, N,N)
#' for (i in 1:(N-1)){
#'   for (j in (i+1):N){
#'      sim_matr[i,j] <- sum(a[,i]==a[,j])/H
#'      sim_matr[j,i] <- sim_matr[i,j]
#'      }
#'}
#'
#' # Obtain a clustering solution via kmeans with multiple random seeds
#'
#' cl <- KMeans(x, centers)$cluster
#'
#' # Find three pivots
#'
#' mus_alg <- MUS(C = sim_matr, clusters = cl)
#'}
#'
#'@export


#############################################
#MUS algorithm
###########################################

MUS <- function(C, clusters, prec_par=10) {

  tol=0
  if (missing(prec_par)) {
    prec_par <- min(10, min(table(clusters)))
  }else {
    prec_par <- min(prec_par, min(table(clusters)))
  }

  nCl=length(unique(clusters))


  if(nCl>=5){stop("The maximum number of clusters allowed is four!")}


  coppie <- which(C <= tol, arr.ind = TRUE)

  ordine <- sort(table(coppie), decreasing = TRUE)

  lista_ord <- as.double(names(ordine))

  subgroups <- matrix(NA, nCl, prec_par)

  for (g in 1:length(unique(clusters))) {

    if (sum(clusters==g)==1) {
      warning("Singleton cluster given as input!")
      geterrmessage()
    }

    subgroups[g, ] <- lista_ord[clusters[lista_ord] == g][1:prec_par]

    if (all(is.na(subgroups[g, ])) == TRUE) {
      stop("The set of candidate pivots is empty, try increasing prec_par!")
    }

  }

  P <- c()
  insieme_card_p <- c()

  p_star <- as.vector(t(subgroups))
  p_star <- p_star[is.na(p_star)==FALSE]
  contatore <- rep(0, length(p_star))

  custom <- function(m, Coa, Pset, clust) {
    i <- 0
    rip_matr <- matrix(NA, 1, 3)
    lista_rip <- list()

    for (p in 1:length(Pset)) {
      if (identical(Coa[c(m, Pset[p]), c(m, Pset[p])], diag(3)) && clust[Pset[p]] != clust[m]) {
        i <- i + 1
        lista_rip[[i]] <- sort(c(m, Pset[p]))
      }
    }
    if (i != 0) {
      rip_matr <- matrix(NA, i, 3)
      for (j in 1:i) {
        rip_matr[j, ] <- lista_rip[[j]]
      }
    }
    return(rip_matr = rip_matr)
  }


  #  clusters[p_star]

  if (nCl == 2) {
    for (u in 1 : length(p_star)) {
      P <- coppie[, 2][coppie[, 1] == p_star[u] & clusters[coppie[, 2]] != clusters[p_star[u]]]
      card_p <- length(P)
      insieme_card_p[u] <- card_p

      if (insieme_card_p[u] > 1) {
        for (p in 1:(insieme_card_p[u] - 1)) {
          if (C[P[p], p_star[u]] <= tol && clusters[P[p]] != clusters[p_star[u]]) {
            contatore[u] = contatore[u] + 1
          }
        }
      }
    }
  }else if (nCl == 3) {
    for (u in 1: length(p_star)){
      P <- coppie[, 2][coppie[, 1] == p_star[u] & clusters[coppie[, 2]] != clusters[p_star[u]]]
      card_p <- length(P)
      insieme_card_p[u] <- card_p
      if (insieme_card_p[u] == 2) {
        if (C[P[1], P[2]] <= tol && clusters[P[1]] != clusters[P[2]] &&
            clusters[p_star[u]] != clusters[P[2]]) {
          contatore[u] = contatore[u] + 1
        }
      }else if (insieme_card_p[u] > 2) {

        for (p in 1:(insieme_card_p[u] - 1)) {
          for (o in (p + 1):insieme_card_p[u]) {
            if (C[P[p], P[o]] <= tol && clusters[P[p]] !=  clusters[P[o]] && clusters[p_star[u]] !=  clusters[P[o]]) {
              contatore[u] = contatore[u] + 1
            }
          }
        }
      }
    }
  } else if (nCl == 4){

    for (u in 1:length(p_star)) {

      P=p_star[which(clusters[p_star]!=clusters[p_star[u]])][C[p_star[u],
                                                               p_star[which(clusters[p_star]!=clusters[p_star[u]])]]==0]
      insieme_card_p[u] <- length(P)

      if (insieme_card_p[u] > 2) {
        a <- combn(P, 2)[, apply(combn(P, 2), 2, function(x)
          sum((identical(C[x, x], diag(2))) & (length(unique(clusters[x])) == 2))) == 1]

        aa <- t(a)

        if (dim(aa)[1] != 0) {
          fg_matrix <- matrix()
          list_result <- lapply(split(aa,seq(NROW(aa))),custom,  C, P, clusters)
          fg_matrix <- do.call(rbind,list_result)

          fg_matrix <- subset(fg_matrix, apply(fg_matrix, 1, function(x) identical(is.na(x),
                                                                                   rep(FALSE, 3))))
          contatore[u] <- nrow(unique.matrix(fg_matrix, incomparables = FALSE, 1))
        }
      }
    }
  }

  if (sum(contatore)==0) {warning("The submatrix of pivots is not identical!")}

  tabella <- cbind(p_star, cl=clusters[p_star], insieme_card_p, contatore)
  maxima <- rep(NA, length(unique(clusters)))

  for (g in 1:nCl) {
    s=subset(tabella[, 4], tabella[, 2] == g)
    maxima[g] <- subset(tabella, tabella[, 2] == g)[which.max(s), 1]

  }
  return(list(pivots = maxima,
              prec_par = prec_par))
}
