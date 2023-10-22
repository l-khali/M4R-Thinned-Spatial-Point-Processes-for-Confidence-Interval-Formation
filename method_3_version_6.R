source("poisson_simulation_method_1.R")
library(spatstat)

subsets_method6 <- function(data, N, alpha, R=99) {
  
  increment <- 1/sqrt(N)
  a <- 1/N
  npoints <- nrow(as.data.frame(data))
  # using same values for radius as paper for now
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
    
  k_vals_current <- data.frame(rep(c(0), each=15))
  
  for (sin in 1:R) {
    for (n in 1:N) {
      xstart <- runif(1,0,1)
      ystart <- runif(1,0,1)
      xend <- xstart + increment
      yend <- ystart + increment
      
      # using modulus so that points that "wrap around" are next to each other
      if (xend <= 1 & yend <= 1) {
        subregion <- as.data.frame(subset(data, xstart <= x & x < xend & ystart <= y & y < yend))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
      } else if (xend > 1 & yend <= 1) {
        subregion <- as.data.frame(subset(data, ((xstart <= x | x < (xend %% 1)) & (ystart <= y & y < yend))))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart)))
      } else if (xend <= 1 & yend > 1) {
        subregion <- as.data.frame(subset(data, (xstart <= x & x < xend) & ((ystart <= y | y < (yend %% 1)))))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart)))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
      } else {
        subregion <- as.data.frame(subset(data, (xstart <= x | x < (xend %% 1)) & (ystart <= y | y < (yend %% 1))))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
      }
          
      ni <- nrow(subregion)
      window <- owin(c(0, increment), c(0,increment))
      subregion <- as.ppp(subregion, window)
          
      # removing coefficient of K to obtain just the sum
      k <- Kest(subregion, r = seq(0.0, 0.14, 0.01), correction=c("isotropic"))
      k_vals_current <- cbind(k_vals_current, as.data.frame(k)["iso"] * ni * (ni - 1) / a) # don't need the rep here?
    }
  
    k_vals <- cbind(k_vals, rowSums(k_vals_current) / (npoints * (npoints - 1)))
  }
  K_est <- as.data.frame(Kest(data, r = seq(0.0, 0.14, 0.01), correction=c("isotropic")))["iso"]
  lower_approx <- c()
  upper_approx <- c()
  t <- qt(1-alpha/2, N-1)
  
  for (radius in 1:15) {
    var <- var(as.numeric(k_vals[radius,-1]))/N #extra division by N as specified in Diggle 2003 pg 52
    lower_approx <- c(lower_approx, K_est[radius,] - t*sqrt(var/N))
    upper_approx <- c(upper_approx, K_est[radius,] + t*sqrt(var/N))
  }
  
  return(cbind(lower_approx, upper_approx))
}

subsets4v6 <- poisson_simulation_subsets6(50,250,4,0.05)
plot(subsets4v6, ylim=c(0,1), type="l", lty=3)
