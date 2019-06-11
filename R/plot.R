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
#' }
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
                     type = c("chains", "hist") ){
  colori <- c("red", "green", "violet", "blue")

  ### checks

  # par
  list_par <- c("mean", "sd", "weight", "all")
  if (sum(par!=list_par)==4){
    stop(paste("object ", "'", par,"'", " not found.
    Please select one among the following parameters:
    mean, sd, weight, all", sep=""))
  }

  # type
  list_type <- c("chains", "hist")
  if (sum(type!=list_type)==2){
    stop(paste("object ", "'", type,"'", " not found.
    Please select one among the following types:
    chains, hist", sep=""))
  }



  ###

  if (missing(par)){
    par <- "all"
  }

  if (is.vector(y)){
    est <- apply(rel_est$rel_mean,2,median)
  }else{
    est <- t(apply(rel_est$rel_mean, c(2,3), median))
  }

  if (par=="mean"){
    rel <- rel_est$rel_mean
    raw <- mcmc$mcmc_mean_raw
  }else if(par=="sd"){
    rel <- rel_est$rel_sd
    raw <- mcmc$mcmc_sd_raw
  }else if(par=="weight"){
    rel <- rel_est$rel_weight
    raw <- mcmc$mcmc_weight_raw
  }else{
    if (is.vector(y)){
      rel <- array(0, dim =c(3, dim(rel_est$rel_mean)[1],
                             dim(rel_est$rel_mean)[2]))
      raw <- array(0, dim =c(3, dim(mcmc$mcmc_mean_raw)[1],
                             dim(mcmc$mcmc_mean_raw)[2]))
      rel[1,,] <-  rel_est$rel_mean
      rel[2,,] <-  rel_est$rel_sd
      rel[3,,] <-  rel_est$rel_weight
      raw[1,,] <-  mcmc$mcmc_mean_raw
      raw[2,,] <-  mcmc$mcmc_sd_raw
      raw[3,,] <-  mcmc$mcmc_weight_raw
    }else{
      rel <- array(0, dim =c(3, dim(rel_est$rel_mean)[1],
                             dim(rel_est$rel_mean)[3]))
      raw <- array(0, dim =c(3, dim(mcmc$mcmc_mean_raw)[1],
                             dim(mcmc$mcmc_mean_raw)[3]))
      rel[1,,] <-  rel_est$rel_mean[,1,]
      rel[2,,] <-  rel_est$rel_mean[,2,]
      rel[3,,] <-  rel_est$rel_weight
      raw[1,,] <-  mcmc$mcmc_mean_raw[,1,]
      raw[2,,] <-  mcmc$mcmc_mean_raw[,2,]
      raw[3,,] <-  mcmc$mcmc_weight_raw

    }

  }
  n.iter <- rel_est$final_it
  true.means <- mcmc$Mu

  if (type=="chains" ){
    mains <- c("mean", "sd", "weight")
    ylabs <- c(expression(mu), expression(tau), expression(pi))
    if (is.vector(y)){
      if (par=="all"){
        par(mfrow=c(2,3), oma=c(0,0,0,0), mar =c(5,3,2,1))
        k <- dim(raw)[3]

        for (j in 1:3){
          #plot
          matplot(raw[j,,], type="l", xlab="",
                  ylab=ylabs[j], main= paste("Raw ", mains[j],"s", sep=""), cex.lab =1.8,
                  cex.main=1.8)
        }
        for (j in 1:3){
          #plot the relabeled
          matplot(rel[j,,], type="l", xlab="Iterations",
                  ylab=ylabs[j], main= paste("Rel ", mains[j],"s", sep=""),
                  cex.main=1.8, cex.lab =1.8)
        }
        cat("Description: traceplots of the raw MCMC chains and the relabelled chains for all the model parameters: means, sds and weights. Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sampler is not able to distinguish between the components.")

      }else{

        k <- dim(raw)[2]
        par(mfrow=c(1,2), oma=c(0,0,0,0), mar =c(5,3,3,1))
        #plot
        matplot(raw, type="l", xlab="Iterations",
                ylab="", main= paste("Raw ", par,"s", sep=""),
                cex.main=1.8,cex.lab =1.8)
        #plot the relabeled
        matplot(rel, type="l",
                xlab="Iterations",
                ylab="",
                main= paste("Rel ", par,"s", sep=""), cex.main=1.8,
                cex.lab =1.8)

        cat(paste("Description: traceplot of the raw MCMC chains and the relabelled chains for the "), par,"s parameters. Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.", sep="")
      }
    }else{
      if (par=="all"){
        par(mfrow=c(2,2), oma=c(0,0,0,0), mar =c(5,4.6,2,1))
        k <- dim(raw)[3]
         mains = c("mean 1st coord", "mean 2nd coord", "weight")
         ylabs = c(expression(mu[,1]), expression(mu[,2]), expression(pi))
        #
        # for (j in 1:3){
        #   #plot
        #   matplot(raw[j,,], type="l", xlab="",
        #           ylab=ylabs[j], main= paste("Raw ", mains[j], sep=""),
        #           cex.main=1.8, cex.lab =1.8)
        # }
        # for (j in 1:3){
        #   #plot the relabeled
        #   matplot(rel[j,,], type="l", xlab="Iterations",
        #           ylab=ylabs[j], main= paste("Rel ", mains[j], sep=""),
        #           cex.lab =1.8, cex.main=1.8)
        # }

        h=1
        plot(raw[1,,h], raw[2, ,h], col=h, pch =1, bg =h,
             cex.main =1.8, main ="Raw means", cex.lab=1.8, xlab=
               expression(mu[1]), ylab =expression(mu[2]),
             xlim= c(min(raw[1,,]-10), max(raw[1,,])+10),
             ylim= c(min(raw[2,,]-10), max(raw[2,,])+10))
        for (h in 2:k){
          points(raw[1,,h], raw[2, ,h], col=h, pch =1, bg=h)
        }

        matplot(raw[3,,], type="l",
                          ylab=ylabs[3], main= paste("Raw ", mains[3],"s", sep=""),
                           cex.main=1.8, cex.lab =1.8, xlab ="Iterations")

        h=1
        plot(rel[1,,h], rel[2, ,h], col=h, pch =1, bg =h,
             cex.main =1.8, main ="Rel means", cex.lab=1.8, xlab=
               expression(mu[1]), ylab =expression(mu[2]),
             xlim= c(min(rel[1,,]-10), max(rel[1,,])+10),
             ylim= c(min(rel[2,,]-10), max(rel[2,,])+10))
        for (h in 2:k){
          points(rel[1,,h], rel[2, ,h], col=h, pch =1, bg=h)
        }

        matplot(rel[3,,], type="l",
                ylab=ylabs[3], main= paste("Rel ", mains[3],"s", sep=""),
                cex.main=1.8, cex.lab =1.8, xlab ="Iterations")




        cat("Description: traceplots of the raw MCMC chains and the relabelled chains for the model parameters means and weights. Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.")


      }else{
        if (par=="mean"){

          k <- dim(raw)[3]
          par(mfrow=c(1,2), oma=c(0,0,0,0), mar =c(5,4.6,2,1))
          # matplot(raw[,1,], type="l", xlab="",
          #         ylab=expression(mu[1]), main= paste("Raw ", par,"s", sep=""),
          #         cex.lab =1.8, cex.main=1.8 )
          # #plot the second component
          # matplot(raw[,2,], type="l", xlab="",
          #         ylab=expression(mu[2]), main= paste("Raw ", par,"s", sep=""),
          #         cex.lab =1.8, cex.main=1.8)
          #
          # #plot the first relabelled component
          # matplot(rel[,1,],type="l", xlab="Iterations",
          #         ylab=expression(mu[1]),
          #         main= paste("Rel ", par,"s", sep=""),cex.lab =1.8, cex.main=1.8)
          #
          # #plot the second relabelled component
          # matplot(rel[,2,],type="l", xlab="Iterations",
          #         ylab=expression(mu[2]),
          #         main= paste("Rel ", par,"s", sep=""),cex.lab =1.8, cex.main=1.8)

          h=1
          plot(raw[,1,h], raw[, 2,h], col=h, pch =1, bg =h,
               cex.main =1.8, main ="Raw means", cex.lab=1.8, xlab=
                 expression(mu[1]), ylab =expression(mu[2]),
               xlim= c(min(raw[,1,]-10), max(raw[,1,])+10),
               ylim= c(min(raw[,2,]-10), max(raw[,2,])+10))
          for (h in 2:k){
            points(raw[,1,h], raw[,2 ,h], col=h, pch =1, bg=h)
          }

          h=1
          plot(rel[,1,h], rel[,2,h], col=h, pch =1, bg =h,
               cex.main =1.8, main ="Rel means", cex.lab=1.8, xlab=
                 expression(mu[1]), ylab =expression(mu[2]),
               xlim= c(min(rel[,1,]-10), max(rel[,1,])+10),
               ylim= c(min(rel[,2,]-10), max(rel[,2,])+10))
          for (h in 2:k){
            points(rel[,1,h], rel[,2 ,h], col=h, pch =1, bg=h)
          }



          cat(paste("Description: traceplot of the raw MCMC chains and the relabelled chains for the "), par,"s parameters (coordinate 1 and 2). Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.", sep="")



        }else if(par=="weight"){
          par(mfrow=c(1,2), oma=c(0,0,0,0), mar =c(5,4.6,2,1))

          matplot(raw, type="l", xlab="Iterations",
                  ylab=expression(pi), main= paste("Raw ", par,"s", sep=""),
                  cex.lab =1.8, cex.main=1.8)


          matplot(rel,type="l", xlab="Iterations",
                  ylab=expression(pi),
                  main= paste("Rel ", par,"s", sep=""),
                  cex.lab =1.8, cex.main=1.8)


          cat(paste("Description: traceplot of the raw MCMC chains and the relabelled chains for the "), par,"s parameters. Each colored chain corresponds to one of the k distinct parameters of the mixture model. Overlapping chains may reveal that the MCMC sample is not able to distinguish between the components.", sep="")

        }else if (par=="sd"){
          return(cat("No sds available in two dimensions: they are the same for each group"))
        }

      }
    }

  }else if(type=="hist"){
    if (is.vector(y)){

      par(mfrow=c(1,2), mar=c(5,5,4,1))
      hist(y, breaks=40, prob = TRUE,
           main = paste("Raw means"), cex.main =1.8,
           col="navajowhite1", border="navajowhite1", cex.lab =1.8)
      points(colMeans(mcmc$mcmc_mean_raw), rep(0, length(true.means)),
             col="red", pch=21,  bg="red")
      lines(density(y), lty=1, lwd=3, col="blue")
      hist(y, breaks=40, prob = TRUE,
           main= paste("Rel means"), cex.main=1.8,
           col="navajowhite1", border="navajowhite1", cex.lab =1.8)
      points(est, rep(0, length(true.means)),
             col="red", pch=21,  bg="red")
      lines(density(y), lty=1, lwd=3, col="blue")

      cat("Description: histograms of the data along with the estimated posterior means (red points) from raw MCMC and relabelling algorithm. The blue line is the estimated density curve.")

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
      #cols <- (freq2D[-1,-1] + freq2D[-1,-(nbins-1)] + freq2D[-(nbins-1),-(nbins-1)] + freq2D[-(nbins-1),-1])/4
      res <- persp(x.bin, y.bin,
                   freq2D,   xlab="\n\n\nx",
                   ylab="\n\n\ny", zlab="\n\n\nf(x,y)",
                   theta=30, phi=30,
                   expand=0.5, ltheta=1,
                   lphi=1,
                   col = "white",
                     #"navajowhite1",
                   shade = 0.1, ticktype = "detailed",
                   main= paste("Rel means"), cex.main=1.8,
                   cex.lab =1.8)

       points(trans3d(est[,1],
                      est[,2], -2,
                      pmat = res), col = "red", pch = 17,
              cex=1.5)
       # points(trans3d(est[,1],
       #                est[,2], max(freq[,3])+1,
       #                pmat = res), col = "red", pch = 16,
       #        cex=1.5)
       #lines(trans3d(est[,1], est[,2],
        #             max(freq[,3])+1, pmat =res), col="red")

      cat("Description: 3d histogram of the data along with the posterior estimates of the relabelled means (red points)")

    }

  }
  # else if (type=="iter"){
  #   par(mfrow=c(1,1))
  #   barplot(n.iter, main="Proportion of valid iterations", ylim=c(0,1),
  #           xlab="Pivotal criterion", ylab="Prop.", names.arg=c(1:7))
  # }

}










#else if (type=="estimates"){
#   if (length(dim(raw))==2){
#     switch.means <- colMeans(raw)
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
#     points(t(apply(raw, c(2,3), mean)), col="red")
#     #plot relabelled estimates
#     plot(true.means, xlim=c( min(true.means, est)-1,
#       max(true.means, est)+1  ), ylim=c(u1,u2),
#       xlab=expression(mu[1]), ylab=expression(mu[2]),
#       main="Relabelled",  pch=3, bg=2,
#       cex.main = 0.7)
#     points(est, col="red")
#   }
# }
