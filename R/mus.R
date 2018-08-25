#'MUS algorithm
#'
#' Finding the pivotal units through a sequential search in the symmetric matrix C
#'
#' @param C Square symmetrix matrix with value bounded in \code{[0,1]}. For instance, a co-association matrix resulting from clustering ensembles.
#' @param clusters An initial group assignment for the \code{N} statistical units in \code{k} groups.
#' @param prec_par  A precision parameter for exploring a greater number of algorithm solutions. Default value is 5.
#' @details
#' See the vignette.
#'
#' @return
#'
#' \item{\code{maxima}}{ The \code{k} maxima units}
#'
#' @example
#'N <- 20
#'H <- 1000
#'a <- matrix(NA, H, N)
#'
#'for (h in 1:H){
#'    a[h,] <- kmeans(x,centers)$cluster
#'}
#' build the similarity matrix
#' sim_matr <- matrix(1, n,n)
#' for (i in 1:(n-1)){
#'   for (j in (i+1):n){
#'      sim_matr[i,j] <- sum(a[,i]==a[,j])/H
#'      sim_matr[j,i] <- sim_matr[i,j]
#'      }
#'}
#'
#' cl <- KMeans(x, centers)$cluster
#' mus_alg <- MUS(C = sim_matr, clusters = cl, prec_par = 5)
#'
#'
#'


#############################################
#MUS algorithm
###########################################

MUS <- function(C, clusters, prec_par){

  #C: C di similarit?
  #clusters: gruppi partizione clustering
  #prec_par: parametro di precisione. Indica la numerosit? del sottoinsieme da cui
  #pescare i pivots

  #clusters<-c()
  #clusters<-as.vector(clusters_list[[ind_id]])
  coppie<-which(C==0, arr.ind = TRUE) #estraggo le coppie
  ordine<-sort(table(which(C==0, arr.ind = TRUE)), decreasing=TRUE) #ordino i pptenziali pivot
  new_lista_ord<-c()
  lista_ord<-as.double(names(ordine))
  subgroups<-matrix(NA,length(unique(clusters)), prec_par)
  for (g in 1:length(unique(clusters))) {
    subgroups[g,]<-lista_ord[clusters[lista_ord]==g][1:prec_par]
    if (all(is.na(subgroups[g,]))==TRUE){
      subgroups[g, 1:min(prec_par,sum(clusters==g))]<-which(clusters==g)[1:min(prec_par,sum(clusters==g))]
    }
  }
  p_star<-c()
  P<-c()
  insieme_card_p<-c()

  if (length(lista_ord)<= prec_par*length(unique(clusters))){
    tol <-seq(0.01, 0.1,length.out= 10)
    sel_tol<-c()
    for (i in 1:10){
    sel_tol[i]<- sum(C<=tol[i])
    }
    coppie<-which(C==0 | C<=sel_tol[1], arr.ind = TRUE) #estraggo le coppie
    ordine<-sort(table(which(C==0 |C<=sel_tol[1], arr.ind = TRUE)), decreasing=TRUE) #ordino i pptenziali pivot
    new_lista_ord<-c()
    lista_ord<-as.double(names(ordine))
    subgroups<-matrix(NA,length(unique(clusters)), prec_par)
    for (g in 1:length(unique(clusters))) {
      subgroups[g,]<-lista_ord[clusters[lista_ord]==g][1:prec_par]
      if (all(is.na(subgroups[g,]))==TRUE){
        subgroups[g, 1:min(prec_par,sum(clusters==g))]<-which(clusters==g)[1:min(prec_par,sum(clusters==g))]
      }
    }
  }


  contatore<-rep(0,min(length(lista_ord), prec_par*length(unique(clusters))))

  if (length(unique(clusters))==2){
    for (u in 1:min(length(lista_ord), (prec_par*length(unique(clusters))))){
      p_star[u]<-lista_ord[u]
      P<-coppie[,2][coppie[,1]==p_star[u]& clusters[coppie[,2]]!=clusters[p_star[u]]]
      card_p<-length(P)
      insieme_card_p[u]<-card_p
      #step.two<-coppie[,2][coppie[,1]==p_star[u]]
      if (insieme_card_p[u]>1){
        for ( p in 1 :(insieme_card_p[u]-1)){

            if (C[P[p],p_star[u]]==0 && clusters[P[p]]!=clusters[p_star[u]] ){

              contatore[u]=contatore[u]+1
            }

        }

      }
      #u=u+1



    }
  }else if (length(unique(clusters))==3){


    # for (u in 1:min(length(lista_ord),(prec_par*length(unique(clusters)))) ){
    #   p_star[u]<-lista_ord[u]
    #   P<-coppie[,2][coppie[,1]==p_star[u]& clusters[coppie[,2]]!=clusters[p_star[u]]]
    #   card_p<-length(P)
    #   insieme_card_p[u]<-card_p
    #   #step.two<-coppie[,2][coppie[,1]==p_star[u]]if (u<=(prec_par*length(unique(clusters)))  ){
    #
    #   #if (u<=card_p){
    #     if (insieme_card_p[u]==2){
    #
    #       if (C[P[1],P[2]]==0 && clusters[P[1]]!=clusters[P[2]] &&
    #           clusters[p_star[u]]!= clusters[P[2]]){
    #
    #         contatore[u]=contatore[u]+1
    #       }
    #
    #     }else if(insieme_card_p[u]>2){
    #     for ( p in 1 :(insieme_card_p[u]-1)){
    #       for (o in (p+1):insieme_card_p[u]){
    #
    #
    #           if (C[P[p],P[o]]==0 && clusters[P[p]]!=clusters[P[o]] &&
    #               clusters[p_star[u]]!= clusters[P[o]]){
    #
    #           contatore[u]=contatore[u]+1
    #         }
    #
    #       }
    #     }
    #
    #   }
    #   #u=u+1
    # #}else{ break }
    # }

    matrice_ord <- cbind(lista_ord, clusters[lista_ord])
    p_star[1:prec_par] <- matrice_ord[ clusters[lista_ord]== unique(clusters)[1], 1][1:prec_par]
    p_star[(prec_par+1):(2*prec_par)  ] <- matrice_ord[ clusters[lista_ord]== unique(clusters)[2],1] [1:prec_par]
    p_star[((2*prec_par)+1):(3*prec_par)  ] <- matrice_ord[clusters[lista_ord]== unique(clusters)[3],1][1:prec_par]

    for (u in 1:min(length(lista_ord),(prec_par*length(unique(clusters)))) ){


      P<-coppie[,2][coppie[,1]==p_star[u]& clusters[coppie[,2]]!=clusters[p_star[u]]]
      card_p<-length(P)
      insieme_card_p[u]<-card_p
      #step.two<-coppie[,2][coppie[,1]==p_star[u]]if (u<=(prec_par*length(unique(clusters)))  ){

      #if (u<=card_p){
      if (insieme_card_p[u]==2){

        if (C[P[1],P[2]]==0 && clusters[P[1]]!=clusters[P[2]] &&
            clusters[p_star[u]]!= clusters[P[2]]){

          contatore[u]=contatore[u]+1
        }

      }else if(insieme_card_p[u]>2){
        for ( p in 1 :(insieme_card_p[u]-1)){
          for (o in (p+1):insieme_card_p[u]){


            if (C[P[p],P[o]]==0 &&
                clusters[P[p]]!=clusters[P[o]] &&
                clusters[p_star[u]]!= clusters[P[o]]){

              contatore[u]=contatore[u]+1
            }

          }
        }

      }
      #u=u+1
      #}else{ break }
    }


  }else{

    new_lista_ord<-c(subgroups[1,],subgroups[2,], subgroups[3,], subgroups[4,])
    for (u in 1:min(length(lista_ord), (prec_par*length(unique(clusters)))) ){

      #p_star[u]<-new_lista_ord[u]
      p_star<-new_lista_ord[is.na(new_lista_ord)==FALSE]
      P<-coppie[,2][coppie[,1]==p_star[u]& clusters[coppie[,2]]!=clusters[p_star[u]]]
      card_p<-length(P)
      insieme_card_p[u]<-card_p


      if (insieme_card_p[u]>2){

        a<-combn(P,2)[,apply(combn(P,2),2,function(x) sum((identical(C[x, x], diag(2))) & (length(unique(clusters[x]))==2)))==1]
        aa<-t(a)

        if (dim(aa)[1]==0){
          contatore[u]<-0
        }else{

        custom<-function(x)
        {
          #rip_matr<-c(0,0,0)
          y<-0
          rip_matr<-matrix(NA, 1,3)
          lista_rip<-list()
          for (p in 1:length(P)) {

            if (identical(C[c(x, P[p]),c(x,P[p])], diag(3)) && clusters[P[p]]!=clusters[x]){

              y<-y+1

              if (y!=0){
              lista_rip[[y]]<-sort(c(x, P[p]))}
            }}



          if (y!=0){
            rip_matr<-matrix(NA, y,3)


          for (j in 1:y){
            rip_matr[j,]<-lista_rip[[j]]
          }

          new_y<-nrow(t(apply(rip_matr,1,function(x) unique(x))))
          }else{new_y=0}




          return(list(new_y=new_y, rip_matr=rip_matr ) )

        }

        fg<-list()

        for (row in 1:nrow(aa)){
          fg[[row]]<-custom(aa[row,])$rip_matr
        }

        fg_matrix<-matrix(rep(NA,3), ncol=3)
        #fg_matrix<-matrix(fg[[1]], nrow= nrow(fg[[1]]), ncol=3)

        for (r in 1:(length(fg))){
          fg_matrix<-rbind(fg_matrix, fg[[r]])
        }

        fg_matrix<-subset(fg_matrix, apply(fg_matrix,1, function(x)  identical(is.na(x), rep(FALSE,3)  )  ))



        contatore[u]<-nrow(unique.matrix(fg_matrix, incomparables = FALSE, 1))
          #sum(apply(aa,1, custom))
        ###########################

        #sol 2: provata. Lenta! ##################
        # contatore[u]<-sum(apply(combn(P,3),2, function(x) sum(identical(C[x, x], diag(3))
        # & length(unique(clusters[x]))==3)))

        #################################

      } }
      #u=u+1

    }
  }

  tabella<-cbind(p_star, clusters[p_star], insieme_card_p, contatore)



  maxima<-c()

  if (length(unique(tabella[,2]))<length(unique(clusters))){
    maxima<-rep(NA,length(unique(clusters)) )

  return(list(maxima=maxima))
  }else{



  p<-c()
  for (g in 1:length(unique(clusters))){
    maxima[g]<- subset(tabella, tabella[,2]==g)[which.max(subset(tabella[,4], tabella[,2]==g)),1]
  }
  return(list(tabella=tabella,
              maxima=maxima,
              contatore=contatore
              #,
              #new=new_lista_ord[is.na(new_lista_ord)==FALSE
               ))

  }}
