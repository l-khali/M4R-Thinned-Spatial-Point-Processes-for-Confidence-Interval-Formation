library(spatstat)

approximation_method <- function(data, N, alpha) {
  # assume N is a square integer
  # assume data is a ppp object
  
  increment <- 1/sqrt(N)
  # using same values for radius as paper for now
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  variances <- c()
  
  for (xstart in head(seq(0,1,length.out = sqrt(N) + 1),-1)) {
    for (ystart in head(seq(0,1,length.out = sqrt(N) + 1),-1)) {
      subregion <- subset(data, xstart <= x & x < xstart + increment & ystart <= y & y < ystart + increment)
      # isotropic correction used in paper
      k <- Kest(subregion, r = seq(0.0, 0.14, 0.01), correction=c("isotropic"))
      k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
    }
  }
  
  # investigate why some K values ara NaN!!!
  k_vals[is.na(k_vals)] <- 0
  
  # calculation sample variance for each radius
  for (row in 1:nrow(k_vals)) {
    variances <- c(variances, var(as.numeric(k_vals[row,1:N+1])))
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



