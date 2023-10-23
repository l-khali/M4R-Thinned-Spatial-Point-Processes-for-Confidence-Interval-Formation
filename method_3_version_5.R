source("poisson_simulation_method_1.R")
library(spatstat)

subsets_method5 <- function(data, N, alpha, R=99) {
  
  increment <- 1/sqrt(N)
  a <- 1/N
  npoints <- nrow(as.data.frame(data))
  # using same values for radius as paper for now
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  window <- owin(c(0, increment), c(0,increment))
  
  for (sim in 1:R) {
    ni_sum <- 0
    k_vals_current <- data.frame(rep(c(0), each=15))
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
        subregion <- as.data.frame(subset(data, (xstart <= x & x < xend) & (ystart <= y | y < (yend %% 1))))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart)))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
      } else {
        subregion <- as.data.frame(subset(data, (xstart <= x | x < (xend %% 1)) & (ystart <= y | y < (yend %% 1))))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
      }
      
      ni <- nrow(subregion)
      ni_sum <- ni_sum + ni
      subregion <- as.ppp(subregion, window)
      
      # removing coefficient of K to obtain just the sum
      k <- Kest(subregion, r = seq(0.0, 0.14, 0.01), correction=c("isotropic"))
      k_vals_current <- cbind(k_vals_current, as.data.frame(k)["iso"] * ni * (ni - 1) / a)
      k_vals_current <- replace(k_vals_current, is.na(k_vals_current), 0)
    }
    k_vals <- cbind(k_vals, rowSums(k_vals_current) / (ni_sum * (ni_sum - 1)))
  }
  
  K_est <- as.data.frame(Kest(data, r = seq(0.0, 0.14, 0.01), correction=c("isotropic")))["iso"]
  lower_approx <- c()
  upper_approx <- c()
  
  for (radius in 1:15) {
    sorted_k_vals <- sort(as.numeric(k_vals[radius,-1]))
    lower <- sorted_k_vals[round((R+1)*(1-alpha/2), digits=0)]
    upper <- sorted_k_vals[round((R+1)*(alpha/2), digits=0)]
    lower_approx <- c(lower_approx, 2*K_est[radius,] - lower)
    upper_approx <- c(upper_approx, 2*K_est[radius,] - upper)
  }
  return(cbind(lower_approx, upper_approx))
}

#subsets4v5 <- poisson_simulation_subsets5(1000,250,4,0.05)

#subsets16v5 <- poisson_simulation_subsets5(1000,250,16,0.05)

#subsets64v5 <- poisson_simulation_subsets5(1000,250,64,0.05)

#png(file="~/Documents/year_4/m4r_spatial_stats/initial_plots/working/poisson_subsets_1000.png",width=600, height=600)
#plot(subsets4v5[-1,], ylim=c(0.5,1), type="l", lty=3, main="Poisson: subsets")
#lines(subsets16v5[-1,], ylim=c(0.5,1), type="l", lty=2)
#lines(subsets64v5[-1,], ylim=c(0.5,1), type="l", lty=1)
#dev.off()
