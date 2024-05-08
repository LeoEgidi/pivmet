context("arguments errors")
test_that("piv_KMeans recognise errors/warnings",{
  N  <- 620
  k  <- 3
  n1 <- 20
  n2 <- 100
  n3 <- 500
  x  <- matrix(NA, N,2)
  truegroup <- c( rep(1,n1), rep(2, n2), rep(3, n3))

  x[1:n1,] <- rmvnorm(n1, c(1,5), sigma=diag(2))
  x[(n1+1):(n1+n2),] <- rmvnorm(n2, c(4,0), sigma=diag(2))
  x[(n1+n2+1):(n1+n2+n3),] <- rmvnorm(n3, c(6,6), sigma=diag(2))

  # k =0
  expect_error(piv_KMeans(x, centers = 0))
  # wrong alg type
  expect_error(piv_KMeans(x, centers = k, alg.type = "diana"))
  # wrong piv.criterion
  expect_error(piv_KMeans(x, centers = k, piv.criterion = "some"))
  # MUS with k > 4
  expect_error(piv_KMeans(x, centers = 5, piv.criterion = "MUS"))

})
