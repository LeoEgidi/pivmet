## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load, warning =FALSE, message = FALSE-------------------------------
library(pivmet)

## ----fish_hist, fig.align ='center'--------------------------------------
data(fish)
y <- fish[,1]
hist(y, breaks=40, prob = TRUE, cex.lab=1.6,
            main ="Fishery data", cex.main =1.7,
            col="navajowhite1", border="navajowhite1")
lines(density(y), lty=1, lwd=3, col="blue")

