context("arguments errors")
test_that("piv_rel recognise errors/warnings",{
  library(bayesmix)
  data(fish)
  y <- fish[,1]
  k <- 5
  nMC <- 5000
  res <- piv_MCMC(y = y, k = k, nMC = nMC)

  # wrong output
  expect_error(piv_rel(mcmc = k ))

})


