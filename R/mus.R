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
#' @author Leonardo Egidi \url{legidi@units.it}, Roberta Pappadà
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
#' sim_matr <- matrix(1, N,N)
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


# function(C, clusters, prec_par){
#
#   #C: C di similarit?
#   #clusters: gruppi partizione clustering
#   #prec_par: parametro di precisione. Indica la numerosit? del sottoinsieme da cui
#   #pescare i pivots
#
#   #clusters<-c()
#   #clusters<-as.vector(clusters_list[[ind_id]])
#
#   if (missing(prec_par)){
#     prec_par <- min(10, min(table(clusters)))
#   }else{
#     prec_par <- min(prec_par, min(table(clusters)))
#   }
#
#   coppie<-which(C==0, arr.ind = TRUE) #estraggo le coppie
#   ordine<-sort(table(which(C==0, arr.ind = TRUE)), decreasing=TRUE) #ordino i pptenziali pivot
#   new_lista_ord<-c()
#   lista_ord<-as.double(names(ordine))
#   subgroups<-matrix(NA,length(unique(clusters)), prec_par)
#   for (g in 1:length(unique(clusters))) {
#     subgroups[g,]<-lista_ord[clusters[lista_ord]==g][1:prec_par]
#     if (all(is.na(subgroups[g,]))==TRUE){
#       subgroups[g, 1:min(prec_par,sum(clusters==g))]<-which(clusters==g)[1:min(prec_par,sum(clusters==g))]
#     }
#   }
#   p_star<-c()
#   P<-c()
#   insieme_card_p<-c()
#
#   if (length(lista_ord)<= prec_par*length(unique(clusters))){
#     tol <-seq(0.01, 0.1,length.out= 10)
#     sel_tol<-c()
#     for (i in 1:10){
#     sel_tol[i]<- sum(C<=tol[i])
#     }
#     coppie<-which(C==0 | C<=sel_tol[1], arr.ind = TRUE) #estraggo le coppie
#     ordine<-sort(table(which(C==0 |C<=sel_tol[1], arr.ind = TRUE)), decreasing=TRUE) #ordino i pptenziali pivot
#     new_lista_ord<-c()
#     lista_ord<-as.double(names(ordine))
#     subgroups<-matrix(NA,length(unique(clusters)), prec_par)
#     for (g in 1:length(unique(clusters))) {
#       subgroups[g,]<-lista_ord[clusters[lista_ord]==g][1:prec_par]
#       if (all(is.na(subgroups[g,]))==TRUE){
#         subgroups[g, 1:min(prec_par,sum(clusters==g))]<-which(clusters==g)[1:min(prec_par,sum(clusters==g))]
#       }
#     }
#   }
#
#
#   contatore<-rep(0,min(length(lista_ord), prec_par*length(unique(clusters))))
#
#   if (length(unique(clusters))==2){
#     for (u in 1:min(length(lista_ord), (prec_par*length(unique(clusters))))){
#       p_star[u]<-lista_ord[u]
#       P<-coppie[,2][coppie[,1]==p_star[u]& clusters[coppie[,2]]!=clusters[p_star[u]]]
#       card_p<-length(P)
#       insieme_card_p[u]<-card_p
#       #step.two<-coppie[,2][coppie[,1]==p_star[u]]
#       if (insieme_card_p[u]>1){
#         for ( p in 1 :(insieme_card_p[u]-1)){
#
#             if (C[P[p],p_star[u]]==0 && clusters[P[p]]!=clusters[p_star[u]] ){
#
#               contatore[u]=contatore[u]+1
#             }
#
#         }
#
#       }
#       #u=u+1
#
#
#
#     }
#   }else if (length(unique(clusters))==3){
#
#
#     # for (u in 1:min(length(lista_ord),(prec_par*length(unique(clusters)))) ){
#     #   p_star[u]<-lista_ord[u]
#     #   P<-coppie[,2][coppie[,1]==p_star[u]& clusters[coppie[,2]]!=clusters[p_star[u]]]
#     #   card_p<-length(P)
#     #   insieme_card_p[u]<-card_p
#     #   #step.two<-coppie[,2][coppie[,1]==p_star[u]]if (u<=(prec_par*length(unique(clusters)))  ){
#     #
#     #   #if (u<=card_p){
#     #     if (insieme_card_p[u]==2){
#     #
#     #       if (C[P[1],P[2]]==0 && clusters[P[1]]!=clusters[P[2]] &&
#     #           clusters[p_star[u]]!= clusters[P[2]]){
#     #
#     #         contatore[u]=contatore[u]+1
#     #       }
#     #
#     #     }else if(insieme_card_p[u]>2){
#     #     for ( p in 1 :(insieme_card_p[u]-1)){
#     #       for (o in (p+1):insieme_card_p[u]){
#     #
#     #
#     #           if (C[P[p],P[o]]==0 && clusters[P[p]]!=clusters[P[o]] &&
#     #               clusters[p_star[u]]!= clusters[P[o]]){
#     #
#     #           contatore[u]=contatore[u]+1
#     #         }
#     #
#     #       }
#     #     }
#     #
#     #   }
#     #   #u=u+1
#     # #}else{ break }
#     # }
#
#     matrice_ord <- cbind(lista_ord, clusters[lista_ord])
#     p_star[1:prec_par] <- matrice_ord[ clusters[lista_ord]== unique(clusters)[1], 1][1:prec_par]
#     p_star[(prec_par+1):(2*prec_par)  ] <- matrice_ord[ clusters[lista_ord]== unique(clusters)[2],1] [1:prec_par]
#     p_star[((2*prec_par)+1):(3*prec_par)  ] <- matrice_ord[clusters[lista_ord]== unique(clusters)[3],1][1:prec_par]
#
#     for (u in 1:min(length(lista_ord),(prec_par*length(unique(clusters)))) ){
#
#
#       P<-coppie[,2][coppie[,1]==p_star[u]& clusters[coppie[,2]]!=clusters[p_star[u]]]
#       card_p<-length(P)
#       insieme_card_p[u]<-card_p
#       #step.two<-coppie[,2][coppie[,1]==p_star[u]]if (u<=(prec_par*length(unique(clusters)))  ){
#
#       #if (u<=card_p){
#       if (insieme_card_p[u]==2){
#
#         if (C[P[1],P[2]]==0 && clusters[P[1]]!=clusters[P[2]] &&
#             clusters[p_star[u]]!= clusters[P[2]]){
#
#           contatore[u]=contatore[u]+1
#         }
#
#       }else if(insieme_card_p[u]>2){
#         for ( p in 1 :(insieme_card_p[u]-1)){
#           for (o in (p+1):insieme_card_p[u]){
#
#
#             if (C[P[p],P[o]]==0 &&
#                 clusters[P[p]]!=clusters[P[o]] &&
#                 clusters[p_star[u]]!= clusters[P[o]]){
#
#               contatore[u]=contatore[u]+1
#             }
#
#           }
#         }
#
#       }
#       #u=u+1
#       #}else{ break }
#     }
#
#
#   }else{
#
#     new_lista_ord<-c(subgroups[1,],subgroups[2,], subgroups[3,], subgroups[4,])
#     for (u in 1:min(length(lista_ord), (prec_par*length(unique(clusters)))) ){
#
#       #p_star[u]<-new_lista_ord[u]
#       p_star<-new_lista_ord[is.na(new_lista_ord)==FALSE]
#       P<-coppie[,2][coppie[,1]==p_star[u]& clusters[coppie[,2]]!=clusters[p_star[u]]]
#       card_p<-length(P)
#       insieme_card_p[u]<-card_p
#
#
#       if (insieme_card_p[u]>2){
#
#         a<-combn(P,2)[,apply(combn(P,2),2,function(x) sum((identical(C[x, x], diag(2))) & (length(unique(clusters[x]))==2)))==1]
#         aa<-t(a)
#
#         if (dim(aa)[1]==0){
#           contatore[u]<-0
#         }else{
#
#         custom<-function(x)
#         {
#           #rip_matr<-c(0,0,0)
#           y<-0
#           rip_matr<-matrix(NA, 1,3)
#           lista_rip<-list()
#           for (p in 1:length(P)) {
#
#             if (identical(C[c(x, P[p]),c(x,P[p])], diag(3)) && clusters[P[p]]!=clusters[x]){
#
#               y<-y+1
#
#               if (y!=0){
#               lista_rip[[y]]<-sort(c(x, P[p]))}
#             }}
#
#
#
#           if (y!=0){
#             rip_matr<-matrix(NA, y,3)
#
#
#           for (j in 1:y){
#             rip_matr[j,]<-lista_rip[[j]]
#           }
#
#           new_y<-nrow(t(apply(rip_matr,1,function(x) unique(x))))
#           }else{new_y=0}
#
#
#
#
#           return(list(new_y=new_y, rip_matr=rip_matr ) )
#
#         }
#
#         fg<-list()
#
#         for (row in 1:nrow(aa)){
#           fg[[row]]<-custom(aa[row,])$rip_matr
#         }
#
#         fg_matrix<-matrix(rep(NA,3), ncol=3)
#         #fg_matrix<-matrix(fg[[1]], nrow= nrow(fg[[1]]), ncol=3)
#
#         for (r in 1:(length(fg))){
#           fg_matrix<-rbind(fg_matrix, fg[[r]])
#         }
#
#         fg_matrix<-subset(fg_matrix, apply(fg_matrix,1, function(x)  identical(is.na(x), rep(FALSE,3)  )  ))
#
#
#
#         contatore[u]<-nrow(unique.matrix(fg_matrix, incomparables = FALSE, 1))
#           #sum(apply(aa,1, custom))
#         ###########################
#
#         #sol 2: provata. Lenta! ##################
#         # contatore[u]<-sum(apply(combn(P,3),2, function(x) sum(identical(C[x, x], diag(3))
#         # & length(unique(clusters[x]))==3)))
#
#         #################################
#
#       } }
#       #u=u+1
#
#     }
#   }
#
#   tabella<-cbind(p_star, clusters[p_star], insieme_card_p, contatore)
#
#
#
#   maxima<-c()
#
#   if (length(unique(tabella[,2]))<length(unique(clusters))){
#     maxima<-rep(NA,length(unique(clusters)) )
#
#     if (is.na(sum(maxima))){
#       return(print("Warning: try increasing the precision parameter prec_par"))
#     }else{
#
#   return(list(pivots=maxima, prec_par = prec_par))
#     }
#   }else{
#
#   p<-c()
#   for (g in 1:length(unique(clusters))){
#     maxima[g]<- subset(tabella, tabella[,2]==g)[which.max(subset(tabella[,4], tabella[,2]==g)),1]
#   }
#
#
#   if (is.na(sum(maxima))){
#     return(print("Warning: try increasing the precision parameter prec_par"))
#   }else{
#   return(list(#tabella=tabella,
#               pivots=maxima,
#               prec_par = prec_par
#               #,contatore=contatore
#               #,
#               #new=new_lista_ord[is.na(new_lista_ord)==FALSE
#                ))
#         }
#   }}


#' @param tol Optional argument. The tolerance value for \code{C} when testing for similarity/co-association.
#' The lowest the tolerance, the highest the requirement for the dissimilarity. Default is 0.
