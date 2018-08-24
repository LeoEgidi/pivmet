#' Plotting outputs from pivotal relabelling
#'
#' Plot and visualize MCMC outputs, posterior relabelled chains and estimates and diagnostics.
#' @param y Data vector or matrix.
#' @param est Pivotal estimates as provided by \code{piv_rel}.
#' @param type Type of plots required. Choose among: \code{"chains"}, \code{"estimates"}, \code{"estimates_hist"}.
#' @param n.iter Number of valid MCMC iterations as provided by \code{piv rel}.
#' @param switch.means  Post-processed chains.
#' @param true.means An estimate for the true means.
#'
#' @examples
#'
#' Fishery data
#'
#' data(fish)
#' y <- fish[,1]
#' N <- length(y)
#' k <- 5
#' nMC <- 5000
#' res <- piv_MCMC(y, k, nMC)
#' rel <- piv_rel(res$mu_switch, res$groupPost, res$clust_sel,
#'                                Mu=res$Mu,
#'                                nMC = nMC)
#' piv_plot(y, res, rel, "chains")
#'
#' piv_plot(y, res, rel, "estimates")
#'
#' piv_plot(y, res, rel, "estimates_hist")
#'
#'
#' @export


piv_plot <- function(y, res, rel, type ){
  colori <- c("red", "green", "violet", "blue")
  est <- rel$mu_rel_median
  chains <-rel$mu_rel_complete
  mu_switch <- res$mu_switch
  n.iter <- rel$Final_It
  true.means <- res$Mu

if (type=="chains" ){
    if (length(dim(mu_switch))==2){

      k <- dim(mu_switch)[2]
      par(mfrow=c(1,2), oma=c(0,0,0,0))
      #plot
      matplot(mu_switch, type="l", xlab="Iterations",
        ylab=expression(mu), main="Raw MCMC chains",
        cex.main=0.8)
      #plot the relabeled
      matplot(chains, type="l",
        xlab="Iterations",
        ylab=expression(mu),
        main=paste("Rel. chains"), cex.main=0.8)
    }else{
      k <- dim(mu_switch)[3]
      par(mfrow=c(2,2), oma=c(0,0,0,0))
      matplot(mu_switch[,1,], type="l", xlab="Iterations",
        ylab=expression(mu[1]), main="Raw MCMC chains",
        cex.main=0.8 )
      #plot the second component
      matplot(mu_switch[,2,], type="l", xlab="Iterations",
        ylab=expression(mu[2]), main="Raw MCMC chains",
        cex.main=0.8)

      #plot the first relabelled component
      matplot(chains[,1,],type="l", xlab="Iterations",
        ylab=expression(mu[1]),
        main=paste("Rel.chains"), cex.main=0.8)

      #plot the second relabelled component
      matplot(chains[,2,],type="l", xlab="Iterations",
        ylab=expression(mu[2]),
        main=paste("Rel. chains"), cex.main=0.8)
    }

  }
  #else if(type=="chains" & all(pivotal.criterion, c(1:6))){
  #   if (length(dim(mu_switch))==2){
  #
  #     k <- dim(mu_switch)[2]
  #     par(mfrow=c(3,2), oma=c(0,0,0,0))
  #
  #     #plot the relabeled
  #     for (j in 1:6){
  #     matplot(chains[[pivotal.criterion[j]]], type="l",
  #       xlab="Iterations",
  #       ylab=expression(mu),
  #       main=paste("Relabelled chains: method", pivotal.criterion[j]))
  #     }
  #
  #   }else{
  #     par(mfrow=c(3,4))
  #
  #     for (j in 1:6){
  #     #plot the first relabelled component
  #     matplot(chains[[pivotal.criterion[j]]][,1,],type="l",
  #       xlab="Iterations",
  #       ylab=expression(mu[1]),
  #       main=paste("Rel. method ", pivotal.criterion[j]))
  #
  #     #plot the second relabelled component
  #     matplot(chains[[pivotal.criterion[j]]][,2,],type="l",
  #       xlab="Iterations",
  #       ylab=expression(mu[2]),
  #       main=paste("Rel. method ", pivotal.criterion[j]) )
  #     }
  #
  #   }
  #
  #
  #}
  else if (type=="estimates"){
    if (length(dim(mu_switch))==2){
      switch.means <- colMeans(mu_switch)
      par(mfrow=c(1,2), oma=c(0,0,0,0), las=1, yaxt="n")
      # raw estimates
      plot( true.means, rep(0.3,length(true.means)),
        axes = FALSE , ylab="",ylim=c(0,1),
        xlim=c( min(true.means,
          est)-2,
          max(true.means,
            est)+2  ),
        main=paste("Raw MCMC estimates" ), cex.main =0.8)
      points(switch.means,
        rep(0, length(true.means)), col="red")
      axis(1)
      axis(1, col = "black", tcl = 0)
      par(yaxt="n")
      axis(2)
      par(yaxt="s")
      axis(2, c(0,0.3), c("Est.", "True"), col = "white", tcl = 0)

      #relabelled estimates
      plot( true.means, rep(0.3,length(true.means)),
        axes = FALSE , ylab="",ylim=c(0,1),
        xlim=c( min(true.means,
          est)-2,
          max(true.means,
            est)+2  ),
        main=paste("Rel. estimates"), cex.main =0.8)
      points(est, rep(0, length(true.means)),
        col="red")
      axis(1)
      axis(1, col = "black", tcl = 0)
      par(yaxt="n")
      axis(2)
      par(yaxt="s")
      axis(2, c(0,0.3), c("Est.", "True"), col = "white", tcl = 0)
    }else{
      par(mfrow=c(1,2), oma =c(0,0,0,0))
      colori<-c("red", "green", "violet", "blue")

      l1<-(3/2)*min(Mu[,1])-max(Mu[,1])/2+5
      l2<-(3/2)*max(Mu[,1])-min(Mu[,1])/2-5
      u1<-(3/2)*min(Mu[,2])-max(Mu[,2])/2
      u2<-(3/2)*max(Mu[,2])-min(Mu[,2])/2

      #plot the raw MCMC estimates
      plot(Mu, xlim=c( min(true.means, est)-2,
        max(true.means,est)+2  ),
        ylim=c(u1,u2), main="Raw MCMC output",
        xlab=expression(mu[1]), ylab=expression(mu[2]))
      points(t(apply(mu_switch, c(2,3), mean)), col="red")
      #for (j in 1:k)
      # points(output_bayes$mu_switch[,,j], col=colori[j])
      #plot relabelled estimates
      plot(Mu, xlim=c( min(true.means, switch.means)-1,
        max(true.means, switch.means)+1  ), ylim=c(u1,u2),
        xlab=expression(mu[1]), ylab=expression(mu[2]),
        main="Relabelled",  pch=3, bg=2)
      points(est, col="red")
    }
  }
  # else if (type=="estimates" & all(pivotal.criterion, c(1:6))){
  #     if (length(dim(mu_switch))==2){
  #       par(mfrow=c(3,2), oma=c(0,0,0,0))
  #       for (j in 1:6){
  #         plot( true.means, rep(0.3,length(true.means)),
  #           axes = FALSE , ylab="",ylim=c(0,1),
  #           xlim=c( min(true.means, est[,pivotal.criterion[j]])-2,
  #             max(true.means, est[,pivotal.criterion[j]])+2  ),
  #           main=paste("Rel. estimates, method", pivotal.criterion[j] ))
  #         points(est[,pivotal.criterion[j]],
  #           rep(0, length(true.means)), col="red")
  #         axis(1)
  #         axis(1, col = "black", tcl = 0)
  #         par(yaxt="n")
  #         axis(2)
  #         par(yaxt="s")
  #         axis(2, c(0,0.3), c("Est.", "True"), col = "white", tcl = 0)
  #       }
  #     }else{
  #     par(mfrow=c(3,2), oma=c(0,0,0,0))
  #     l1<-(3/2)*min(Mu[,1])-max(Mu[,1])/2+5
  #     l2<-(3/2)*max(Mu[,1])-min(Mu[,1])/2-5
  #     u1<-(3/2)*min(Mu[,2])-max(Mu[,2])/2
  #     u2<-(3/2)*max(Mu[,2])-min(Mu[,2])/2
  #
  #
  #      for (j in 1:6){
  #      plot(Mu, xlim=c(l1,l2), ylim=c(u1,u2),
  #        xlab=expression(mu[1]), ylab=expression(mu[2]),
  #        main=paste("Relabelled - method", pivotal.criterion[j], sep=" "),
  #        pch=3, bg=2)
  #      points(est[,,pivotal.criterion[j]], col="red")
  #      }
  #   }
  #
  #   }
  else if(type=="estimates_hist"){
    if (length(dim(mu_switch))==2){
      par(mfrow=c(1,2))
      hist(y, breaks=40, prob = TRUE,
        main ="Raw MCMC estimates", cex.main =0.8)
      points(colMeans(mu_switch), rep(0, length(true.means)),
        col="red", pch=21,  bg="red")
      hist(y, breaks=40, prob = TRUE,
        main=paste("Rel. estimates" ), cex.main=0.8 )
      points(est, rep(0, length(true.means)),
        col="red", pch=21,  bg="red")
    }else{
      par(mfrow=c(1,1), mar=c(3,3,2,1), oma =c(1,1,1,1))
      #   NBiv_mix <- function(x, y, k, rho, mu_x, mu_y,
      #     sigma_1x, sigma_1y, sigma_2x,sigma_2y, p){
      #     a <- (2*pi*sigma_1x*sigma_1y*sqrt(1-rho^2))^(-1)
      #     a2 <- (2*pi*sigma_2x*sigma_2y*sqrt(1-rho^2))^(-1)
      #     distr <- list()
      #     for (j in 1:k){
      #     distr[[j]] <-
      #     p*a*exp(-.5*(1)*(1-rho^2)^(-1)*
      #         ( ( (x-mu_x[j])/sigma_1x )^2+
      #             ((y-mu_y[j])/sigma_1y)^2
      #             -2*rho*(x-mu_x[j])*(y-mu_y[j])/
      #             (sigma_1x*sigma_1y)   ))+
      #             (1-p)*a2*exp(-.5*(1)*(1-rho^2)^(-1)*
      #             ( ( (x-mu_x[j])/sigma_2x )^2+
      #             ((y-mu_y[j])/sigma_2y)^2
      #             -2*rho*(x-mu_x[j])*(y-mu_y[j])/
      #             (sigma_2x*sigma_2y)   ))
      # }
      #
      #     sum_vec <- matrix(NA, k, length(distr[[j]]))
      #     for (j in 1:k){
      #       sum_vec[j,] <- as.vector(distr[[j]])
      #     }
      #
      #     return(apply(sum_vec,2,sum))
      #
      #     }
      #
      #   xx<-seq(min(y[,1])-1, max(y[,1])+1, length.out = min(100, length(y[,1])/2))
      #   yy<-seq(min(y[,2])-1, max(y[,2])+1, length.out = min(100, length(y[,1])/2))
      #   mu_x=add[,1]
      #   mu_y=add[,2]
      #   # poniamo rho=1/2
      #   par(mfrow=c(1,1), oma= c(0,0,0,0))
      #   z<-outer(xx, yy, NBiv_mix, k = length(add[,1]),
      #     mu_x=mu_x, mu_y=mu_y, sigma_1x= add2[1], sigma_1y=add2[2],
      #     sigma_2x=add2[3], sigma_2y=add2[4], p =add2[5], rho=0.5)
      #   #Raw MCMC
      #   persp(xx, yy, z, theta=30, phi=30, xlab="x", ylab="y", zlab="f(x,y)",
      #     expand=0.5, ltheta=120,
      #     col = "lightblue", shade = 0.1, ticktype = "detailed" ) -> res
      #   points(trans3d(t(apply(mu_switch, c(2,3), mean))[,1],
      #     t(apply(mu_switch, c(2,3), mean))[,2], 0, pmat = res), col = 2, pch = 16)
      #   #Relabelled
      #    persp(xx, yy, z, theta=30, phi=30, xlab="x", ylab="y", zlab="f(x,y)",
      #     expand=0.5, ltheta=120,
      #      col = "lightblue",
      #      shade = 0.1, ticktype = "detailed",
      #      main=
      #        paste("Rel. estimates: method ",
      #          pivotal.criterion), cex.main=0.8) -> res
      #    points(trans3d(est[,1,pivotal.criterion],
      #      est[,2,pivotal.criterion], 0, pmat = res), col = 2, pch = 16)

      # 3d histogram

      xy <- y
      nbins <- 20
      x.bin <- seq(floor(min(xy[,1])),
        ceiling(max(xy[,1])), length=nbins)
      y.bin <- seq(floor(min(xy[,2])),
        ceiling(max(xy[,2])), length=nbins)
      freq <-  as.data.frame(table(findInterval(xy[,1],
        x.bin),findInterval(xy[,2], y.bin)))
      freq[,1] <- as.numeric(freq[,1])
      freq[,2] <- as.numeric(freq[,2])
      freq2D <- diag(nbins)*0
      freq2D[cbind(freq[,1], freq[,2])] <- freq[,3]

      #par(mfrow=c(1,2))
      #image(x.bin, y.bin, freq2D,
      #col=topo.colors(max(freq2D)))
      #contour(x.bin, y.bin, freq2D,
      #add=TRUE, col=rgb(1,1,1,.7))
      #palette(rainbow(max(freq2D)))
      #cols <- (freq2D[-1,-1] +
      #freq2D[-1,-(nbins-1)] +
      # freq2D[-(nbins-1),-(nbins-1)] +
      #freq2D[-(nbins-1),-1])/4
      res <- persp(x.bin, y.bin,
        freq2D, theta=30, phi=30, xlab="\n\n\nx",
        ylab="\n\n\ny", zlab="\n\n\nf(x,y)",
        expand=0.5, ltheta=120,
        col = "lightblue",
        shade = 0.1, ticktype = "detailed",
        main=
          paste("Rel. estimates"), cex.main=0.8)
      points(trans3d(est[,1],
        est[,2], 0,
        pmat = res), col = "red", pch = 16)


    }

  }
  # else if(type=="estimates_hist" & all(pivotal.criterion, c(1:6))){
  #     par(mfrow=c(3,2))
  #     hist(y, breaks=40, prob = TRUE,
  #       main=
  #       paste("Real data and relabelled estimates (", pivotal.criterion[1], ")" ) )
  #       points(est[,pivotal.criterion[1]], rep(0, length(true.means)),
  #       col="red", pch=21,  bg="red")
  #
  #     for (j in 2:6){
  #     hist(y, breaks=40, prob = TRUE,
  #     main=
  #     paste("Real data and relabelled estimates (", pivotal.criterion[j], ")" ) )
  #     points(est[,pivotal.criterion[j]], rep(0, length(true.means)),
  #       col="red", pch=21,  bg="red")
  #     }
  #
  #     }
  else if (type=="iter"){
    par(mfrow=c(1,1))
    barplot(n.iter, main="Proportion of valid iterations", ylim=c(0,1),
      xlab="Pivotal criterion", ylab="Prop.", names.arg=c(1:7))
  }

}

#   par(mfrow=c(1,2), oma =c(0,0,0,0))
#   colori<-c("red", "green", "violet", "blue")
#
# l1<-(3/2)*min(Mu[,1])-max(Mu[,1])/2+5
# l2<-(3/2)*max(Mu[,1])-min(Mu[,1])/2-5
# u1<-(3/2)*min(Mu[,2])-max(Mu[,2])/2
# u2<-(3/2)*max(Mu[,2])-min(Mu[,2])/2
#
# plot(Mu, xlim=c(l1,l2), ylim=c(u1,u2), main="Raw MCMC-Pre Processing", xlab=expression(mu[1]), ylab=expression(mu[2]))
# for (j in 1:k)
#   points(output_bayes$mu_pre_switch_compl[,,j], col=colori[j])
#
#
# plot(Mu, xlim=c(l1,l2), ylim=c(u1,u2), main="Raw MCMC output", xlab=expression(mu[1]), ylab=expression(mu[2]))
# for (j in 1:k)
#   points(output_bayes$mu_switch[,,j], col=colori[j])
#
#
#
#
# plot(Mu, xlim=c(l1,l2), ylim=c(u1,u2),xlab=expression(mu[1]),
#   ylab=expression(mu[2]), main="Raw MCMC output", col="black", pch=3, bg=2)
# points(t(apply(output_bayes$mu_switch, c(2,3), mean)), col="red")
#

