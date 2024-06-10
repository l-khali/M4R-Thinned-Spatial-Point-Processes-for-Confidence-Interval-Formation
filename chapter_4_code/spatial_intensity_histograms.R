thinning_param <- 0.8
normalised_intensity <- c()
is <- c()
for (sim in 1:1) {
  p <- rpoispp(500)
  i_hat <- npoints(p)/1
  for (thin_sim in 1:500) {
    subp <- rthin(p,thinning_param)
    is <- c(is, npoints(subp)/thinning_param)
    normalised_intensity <- c(normalised_intensity, (npoints(subp)/thinning_param - i_hat))
  }
}

normalised_intensity2 <- normalised_intensity/sqrt(var(is)*thinning_param)

h1 <- hist(normalised_intensity2)
h2 <- hist(rnorm(500,0,1))
plot(h1,col=3)
plot(h2,col=2, add=TRUE)

