#' Plotting outputs from pivotal relabelling
#'
#' Plot and visualize MCMC outputs, posterior relabelled chains and estimates and diagnostics.
#' @param y Data vector or matrix.
#' @param mcmc The ouptut of the raw MCMC sampling, as provided by \code{piv_MCMC}.
#' @param rel_est Pivotal estimates as provided by \code{piv_rel}.
#' @param type Type of plots required. Choose among: \code{"chains"}, \code{"estimates"}, \code{"hist"}.
#'
#' @examples
#'
#' # Fishery data
#'\dontrun{
#' data(fish)
#' y <- fish[,1]
#' N <- length(y)
#' k <- 5
#' nMC <- 5000
#' res <- piv_MCMC(y = y, k = k, nMC = nMC)
#' rel <- piv_rel(mcmc=res, nMC = nMC)
#' piv_plot(y, res, rel, "chains")
#' piv_plot(y, res, rel, "estimates")
#' piv_plot(y, res, rel, "hist")
#' {}
#'
#' @author
#'
#' Leonardo Egidi \url{legidi@units.it}
#'
#' @export


piv_plot <- function(y,
                     mcmc,
                     rel_est,
                     par = c("mean", "sd", "weight", "all"),
                     type = c("chains", "estimates", "hist") ){
  colori <- c("red", "green", "violet", "blue")

  if (is.vector(y)){
    est <- apply(rel_est$rel_mean,2,median)
  }else{
    est <- apply(rel_est$rel_mean, c(2,3), median)
  }

  if (par=="mean"){
  chains <- rel_est$rel_mean
  mu_switch <- mcmc$mcmc_mean
  }else if(par=="sd"){
  chains <- rel_est$rel_sd
  mu_switch <- mcmc$mcmc_sd
  }else if(par=="weight"){
    chains <- rel_est$rel_weight
    mu_switch <- mcmc$mcmc_weight
  }
  n.iter <- rel_est$final_it
  true.means <- mcmc$Mu

if (type=="chains" ){
    if (length(dim(mu_switch))==2){

      k <- dim(mu_switch)[2]
      par(mfrow=c(1,2), oma=c(0,0,0,0), mar =c(5,4,2,1))
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
      par(mfrow=c(2,2), oma=c(0,0,0,0), mar =c(5,4,2,1))
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
  #else if (type=="estimates"){
  #   if (length(dim(mu_switch))==2){
  #     switch.means <- colMeans(mu_switch)
  #     par(mfrow=c(1,2), oma=c(0,0,0,0),
  #         mar=c(5,4,2,1), las=1, yaxt="n")
  #     # raw estimates
  #     plot( true.means, rep(0.3,length(true.means)),
  #       axes = FALSE , ylab="",ylim=c(0,1),
  #       xlim=c( min(true.means,
  #         est)-2,
  #         max(true.means,
  #           est)+2  ),
  #       main=paste("Raw MCMC estimates" ), cex.main =0.8)
  #     points(switch.means,
  #       rep(0, length(true.means)), col="red")
  #     axis(1)
  #     axis(1, col = "black", tcl = 0)
  #     par(yaxt="n")
  #     axis(2)
  #     par(yaxt="s")
  #     axis(2, c(0,0.3), c("Est.", "True"), col = "white", tcl = 0)
  #
  #     #relabelled estimates
  #     plot( true.means, rep(0.3,length(true.means)),
  #       axes = FALSE , ylab="",ylim=c(0,1),
  #       xlim=c( min(true.means,
  #         est)-2,
  #         max(true.means,
  #           est)+2  ),
  #       main=paste("Rel. estimates"), cex.main =0.8)
  #     points(est, rep(0, length(true.means)),
  #       col="red")
  #     axis(1)
  #     axis(1, col = "black", tcl = 0)
  #     par(yaxt="n")
  #     axis(2)
  #     par(yaxt="s")
  #     axis(2, c(0,0.3), c("Est.", "True"), col = "white", tcl = 0)
  #   }else{
  #     par(mfrow=c(1,2), oma =c(0,0,0,0),  mar=c(5,4,2,0.7))
  #     colori<-c("red", "green", "violet", "blue")
  #
  #     l1<-(3/2)*min(true.means[,1])-max(true.means[,1])/2+5
  #     l2<-(3/2)*max(true.means[,1])-min(true.means[,1])/2-5
  #     u1<-(3/2)*min(true.means[,2])-max(true.means[,2])/2
  #     u2<-(3/2)*max(true.means[,2])-min(true.means[,2])/2
  #
  #     #plot the raw MCMC estimates
  #     plot(true.means, xlim=c( min(true.means, est)-2,
  #       max(true.means,est)+2  ),
  #       ylim=c(u1,u2), main="Raw MCMC",
  #       xlab=expression(mu[1]), ylab=expression(mu[2]), pch =3,
  #       cex.main = 0.7)
  #     points(t(apply(mu_switch, c(2,3), mean)), col="red")
  #     #plot relabelled estimates
  #     plot(true.means, xlim=c( min(true.means, est)-1,
  #       max(true.means, est)+1  ), ylim=c(u1,u2),
  #       xlab=expression(mu[1]), ylab=expression(mu[2]),
  #       main="Relabelled",  pch=3, bg=2,
  #       cex.main = 0.7)
  #     points(est, col="red")
  #   }
  # }
  else if(type=="hist"){
    if (length(dim(mu_switch))==2){
      par(mfrow=c(1,2), mar=c(5,5,4,2))
      hist(y, breaks=40, prob = TRUE,
        main ="Raw MCMC estimates", cex.main =0.8,
        col="navajowhite1", border="navajowhite1")
      points(colMeans(mu_switch), rep(0, length(true.means)),
        col="red", pch=21,  bg="red")
      lines(density(y), lty=1, lwd=3, col="blue")
      hist(y, breaks=40, prob = TRUE,
        main=paste("Rel. estimates" ), cex.main=0.8,
        col="navajowhite1", border="navajowhite1")
      points(est, rep(0, length(true.means)),
        col="red", pch=21,  bg="red")
      lines(density(y), lty=1, lwd=3, col="blue")
    }else{
      par(mfrow=c(1,1), mar=c(3,3,2,1), oma =c(0,0,0,0))
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
      res <- persp(x.bin, y.bin,
        freq2D,   xlab="\n\n\nx",
        ylab="\n\n\ny", zlab="\n\n\nf(x,y)",
        theta=30, phi=30,
        expand=0.5, ltheta=120,
        col = "navajowhite1",
        shade = 0.1, ticktype = "detailed",
        main=
          paste("Rel. estimates"), cex.main=0.8)
      points(trans3d(est[,1],
        est[,2], 0,
        pmat = res), col = "red", pch = 16)


    }

  }else if (type=="iter"){
    par(mfrow=c(1,1))
    barplot(n.iter, main="Proportion of valid iterations", ylim=c(0,1),
      xlab="Pivotal criterion", ylab="Prop.", names.arg=c(1:7))
  }

}

