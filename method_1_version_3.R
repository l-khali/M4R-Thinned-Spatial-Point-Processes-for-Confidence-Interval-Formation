library(spatstat)

approximation_method_3 <- function(data, N, alpha) {
  # assume N is a square integer
  # assume data is a ppp object
  
  increment <- 1/sqrt(N)
  # using same values for radius as paper for now
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  variances <- c()
  
  for (xstart in head(seq(0,1,length.out = sqrt(N) + 1),-1)) {
    for (ystart in head(seq(0,1,length.out = sqrt(N) + 1),-1)) {
      window <- owin(c(xstart,xstart+increment),c(ystart,ystart+increment))
      df <- as.data.frame(subset(data, xstart <= x & x < xstart + increment & ystart <= y & y < ystart + increment))
      if (nrow(df)==0) {
        k_vals <- cbind(k_vals, rep(c(0), each=length(r)))
      } else {
        #isotropic correction used in paper
        k <- Kest(data, r = seq(0.0, 0.14, 0.01), correction=c("isotropic"), domain=window)
        k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
      }
    }
  }
  
  
  # calculation sample variance for each radius
  for (row in 1:nrow(k_vals)) {
    variances <- c(variances, var(as.numeric(k_vals[row,-1])))
  }
  
  # ripley's K of entire region for each radius
  K_est <- as.data.frame(Kest(data, r = seq(0.0, 0.14, 0.01), correction=c("isotropic")))["iso"]
  # t-value
  t <- qt((1-alpha/2), df = N-1)
  
  lower_approx <- c()
  upper_approx <- c()
  for (i in 1:15) {
    lower_approx <- c(lower_approx, K_est[i,] - t * sqrt(variances[i]/N))
    upper_approx <- c(upper_approx, K_est[i,] + t * sqrt(variances[i]/N))
  }
  
  K_actual <- rep(c(pi),each=15) * r * r
  result <- cbind(r, lower_approx, upper_approx)
  return(cbind(lower_approx, upper_approx))
}

coverage_approximation_4 <- poisson_simulation_approximation3(1000,250,4,0.05)
plot(coverage_approximation_4[-1,], ylim = c(0.5,1), main="Poisson: splitting", type="l", lty=3)

coverage_approximation_16 <- poisson_simulation_approximation3(1000,250,16,0.05)
lines(coverage_approximation_16[-1,], ylim = c(0.5,1), main="Approximation, N=16", type="l", lty=2)

coverage_approximation_64 <- poisson_simulation_approximation3(1000,250,64,0.05)
lines(coverage_approximation_64[-1,], ylim = c(0.5,1), main="Approximation, N=64")


