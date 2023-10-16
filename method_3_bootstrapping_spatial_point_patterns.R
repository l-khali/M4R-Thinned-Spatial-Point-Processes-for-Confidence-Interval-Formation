source("poisson_simulation_method_1.R")
library(spatstat)

subsets_method <- function(data, N, alpha, R=99) {
  
  increment <- 1/sqrt(N)
  # using same values for radius as paper for now
  r <- seq(0.0, 0.14, 0.01)
  # K_actual <- rep(c(pi),each=15) * r * r
  k_vals <- data.frame(r)
  
  for (sim in 1:R) {
    
    xstart <- runif(1,0,1)
    ystart <- runif(1,0,1)
    xend <- xstart + 1/sqrt(N)
    yend <- ystart + 1/sqrt(N)
    
    # using modulus so that points that "wrap around" are next to each other
    if (xend <= 1 & yend <= 1) {
      subregion <- as.data.frame(subset(data, xstart <= x & x < xend & ystart <= y & y < yend))
      subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
      subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
    } else if (xend > 1 & yend <= 1) {
      subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & ystart <= y & y < yend))
      subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
      subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
    } else if (xend <= 1 & yend > 1) {
      subregion <- as.data.frame(subset(data, xstart <= x & x < xend & (ystart <= y | y < yend %% 1)))
      subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
      subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
    } else {
      subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y | y < yend %% 1)))
      subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
      subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
    }
    
    window <- owin(c(0, increment), c(0,increment))
    subregion <- as.ppp(subregion, window)
    plot(subregion)
    
    k <- Kest(subregion, r = seq(0.0, 0.14, 0.01), correction=c("isotropic"))
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

# performs better than paper
coverage4 <- poisson_simulation_subsets(1000, 250, 4, 0.05)
coverage16 <- poisson_simulation_subsets(1000, 250, 16, 0.05)
# error when N is very large
coverage64 <- poisson_simulation_subsets(1000, 250, 64, 0.05)

plot(coverage4, ylim=c(0,1), main="Coverage for subsets, N=4")
plot(coverage16, ylim=c(0,1), main="Coverage for subsets, N=16")
plot(coverage64, ylim=c(0,1), main="Coverage for subsets, N=64")
