library(spatstat)
source("poisson_simulation_method_1.R")

approximation_method_2 <- function(data, N, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  v <- varblock(data, Kest, nx=sqrt(N), ny=sqrt(N), confidence=1-alpha, r=r)
  lower <- v$loiso
  upper <- v$hiiso
  print(paste0("lower:",lower))
  print(paste0("upper:",upper))
  return(cbind(lower,upper))
}

#coverage_approximation_4 <- poisson_simulation_approximation2(1000,250,4,0.05)
#plot(coverage_approximation_4, ylim = c(0,1), main="Poisson: splitting", type="l", lty=3)

#coverage_approximation_16 <- poisson_simulation_approximation2(1000,250,16,0.05)
#lines(coverage_approximation_16, ylim = c(0,1), main="Approximation, N=16", type="l", lty=2)

#coverage_approximation_64 <- poisson_simulation_approximation2(1000,250,64,0.05)
#lines(coverage_approximation_64, ylim = c(0,1), main="Approximation, N=64")

