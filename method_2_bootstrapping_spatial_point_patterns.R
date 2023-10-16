source("poisson_simulation_method_1.R")
library(spatstat)

tiling_method <- function(data, N, alpha, R=99) {
  
  increment <- 1/sqrt(N)
  # using same values for radius as paper for now
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  k_vals <- data.frame(r)
  
  # data frame of increments to use for repositioning subregions
  increments <- data.frame()
  for (i in 1:sqrt(N)) {
    for (j in 1:sqrt(N)) {
      increments <- rbind(increments, c((i-1)*increment, (j-1)*increment))
    }
  }
  
  for (sim in 1:R) {
    
    new_process <- data.frame()
    
    for (n in 1:N) {
      
      xstart <- runif(1,0,1)
      ystart <- runif(1,0,1)
      xend <- xstart + 1/sqrt(N)
      yend <- ystart + 1/sqrt(N)
      
      if (xend <= 1 & yend <= 1) {
        subregion <- as.data.frame(subset(data, xstart <= x & x < xend & ystart <= y & y < yend))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[n,1])
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[n,2])
        new_process <- rbind(new_process, subregion)
      } else if (xend > 1 & yend <= 1) {
        subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & ystart <= y & y < yend))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[n,1])
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[n,2])
        new_process <- rbind(new_process, subregion)
      } else if (xend <= 1 & yend > 1) {
        subregion <- as.data.frame(subset(data, xstart <= x & x < xend & (ystart <= y | y < yend %% 1)))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[n,1])
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[n,2])
        new_process <- rbind(new_process, subregion)
      } else {
        subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y | y < yend %% 1)))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[n,1])
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[n,2])
        new_process <- rbind(new_process, subregion)
      }
    }
    
    window <- owin( c(0, 1), c(0,1) )
    new_ppp <- as.ppp(new_process, window)
    
    k <- Kest(new_ppp, r = seq(0.0, 0.14, 0.01), correction=c("isotropic"))
    k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
    
  }
  
  K_est <- as.data.frame(Kest(data, r = seq(0.0, 0.14, 0.01), correction=c("isotropic")))["iso"]
  lower_approx <- c()
  upper_approx <- c()
  
  for (radius in 1:15) {
    sorted_k_vals <- sort(as.numeric(k_vals[radius,-1]))
    lower <- sorted_k_vals[(R+1)*(1-alpha/2)]
    upper <- sorted_k_vals[(R+1)*(alpha/2)]
    lower_approx <- c(lower_approx, 2*K_est[radius,] - lower)
    upper_approx <- c(upper_approx, 2*K_est[radius,] - upper)
  }
  return(cbind(lower_approx, upper_approx))
}

# performing similarly to the paper
coverage <- poisson_simulation_tiling(1000, 250, 4, 0.05)
plot(coverage, ylim=c(0,1), main="Coverage for tiling, N=4")

coverage16 <- poisson_simulation_tiling(1000, 250, 16, 0.05)
plot(coverage16, ylim=c(0,1), main="Coverage for tiling, N=16")

coverage64 <- poisson_simulation_tiling(1000, 250, 64, 0.05)
plot(coverage64, ylim=c(0,1), main="Coverage for tiling, N=64")

#data <- rpoispp(250)
#tiling_results <- tiling_method(data, 4, 0.05, 100)
#r <- seq(0.0, 0.14, 0.01)
#K_actual <- rep(c(pi),each=15) * r * r
#coverage <- coverage_probability(tiling_results, K_actual)
#plot(coverage)

#plot(r, tiling_results[,1])
#lines(r, tiling_results[,2])
#lines(r, K_actual)
  
