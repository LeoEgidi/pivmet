context("arguments errors")
test_that("piv_sel recognise errors/warnings",{
  data(iris)
  # select the columns of variables
  x<- iris[,1:4]
  N <- nrow(x)
  H <- 1000
  a <- matrix(NA, H, N)

  # Perform H k-means partitions

  for (h in 1:H){
    a[h,] <- kmeans(x, centers = 3)$cluster
  }
  # Build the co-association matrix

  C <- matrix(NA, N,N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      C[i,j] <- sum(a[,i]==a[,j])/H
      C[j,i] <- C[i,j]
    }}

  km <- kmeans(x, centers =3)

  data(iris)
# select the columns of variables
x<- iris[,1:4]
N <- nrow(x)
H <- 1000
a <- matrix(NA, H, N)

# Perform H k-means partitions

for (h in 1:H){
 a[h,] <- kmeans(x, centers = 3)$cluster
}
# Build the co-association matrix

C <- matrix(NA, N,N)
for (i in 1:(N-1)){
 for (j in (i+1):N){
   C[i,j] <- sum(a[,i]==a[,j])/H
   C[j,i] <- C[i,j]
 }}

km <- kmeans(x, centers =3)


  # wrong number of cluster indicators
  expect_error(piv_sel(C, clusters = 2))

  # not a co-association matrix
  C[2,1] = 4
  expect_error(piv_sel(C, clusters = km$cluster))

})
