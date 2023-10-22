source("poisson_simulation_method_1.R")
library(spatstat)

subsets_method4 <- function(data, N, alpha, R=99) {
  
  increment <- 1/sqrt(N)
  a <- 1/N
  npoints <- nrow(as.data.frame(data))
  # using same values for radius as paper for now
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  
  for (sim in 1:R) {
    k_vals_current <- data.frame(rep(c(0), each=15))
    for (n in 1:N) {
      xstart <- runif(1,0,1)
      ystart <- runif(1,0,1)
      xend <- xstart + increment
      yend <- ystart + increment
      
      # using modulus so that points that "wrap around" are next to each other
      if (xend <= 1 & yend <= 1) {
        subregion <- as.data.frame(subset(data, xstart <= x & x < xend & ystart <= y & y < yend))
        window <- owin(c(xstart,xend),c(ystart,yend))
      } else if (xend > 1 & yend <= 1) {
        subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & ystart <= y & y < yend))
        window1 <- owin(c(xstart, 1), c(ystart, yend))
        window2 <- owin(c(0,xend%%1), c(ystart, yend))
        window <- union.owin(window1, window2)
      } else if (xend <= 1 & yend > 1) {
        subregion <- as.data.frame(subset(data, xstart <= x & x < xend & (ystart <= y | y < yend %% 1)))
        window1 <- owin(c(xstart, xend), c(ystart, 1))
        window2 <- owin(c(xstart,xend), c(0, yend%%1))
        window <- union.owin(window1, window2)
      } else {
        subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y | y < yend %% 1)))
        window1 <- owin(c(xstart, 1), c(ystart, 1))
        window2 <- owin(c(0,xend%%1), c(ystart, 1))
        window3 <- owin(c(0,xend%%1), c(0,yend%%1))
        window4 <- owin(c(xstart, 1), c(0,yend%%1))
        window <- union.owin(window1,window2,window3,window4)
      }
      
      ni <- nrow(subregion)
      subregion <- as.ppp(subregion, window)
      
      # removing coefficient of K to obtain just the sum
      k <- Kest(subregion, r = seq(0.0, 0.14, 0.01), correction=c("isotropic"))
      k_vals_current <- cbind(k_vals_current, as.data.frame(k)["iso"] * rep(c(ni * (ni - 1) / a),each=15)) # don't need the rep here?
    }
    k_vals <- cbind(k_vals, rowSums(k_vals_current) / (npoints * (npoints - 1)))
  }
  
  K_est <- as.data.frame(Kest(data, r = seq(0.0, 0.14, 0.01), correction=c("isotropic")))["iso"]
  lower_approx <- c()
  upper_approx <- c()
  
  for (radius in 1:15) {
    sorted_k_vals <- sort(as.numeric(k_vals[radius,-1]))
    lower <- sorted_k_vals[round(R+1, digits=0)*(1-alpha/2)]
    upper <- sorted_k_vals[round(R+1, digits=0)*(alpha/2)]
    lower_approx <- c(lower_approx, 2*K_est[radius,] - lower)
    upper_approx <- c(upper_approx, 2*K_est[radius,] - upper)
  }
  return(cbind(lower_approx, upper_approx))
}

coverage_subsets_4_v4 <- poisson_simulation_subsets4(100,250,4,0.05)
plot(coverage_subsets_4_v4, ylim=c(0.5,1), type="l", lty=3)

coverage_subsets_16_v4 <- poisson_simulation_subsets4(100,250,16,0.05)
lines(coverage_subsets_16_v4, ylim=c(0.5,1), type="l", lty=2)

coverage_subsets_4_v4_10 <- poisson_simulation_subsets4(100,50,4,0.05)
plot(coverage_subsets_4_v4_10, ylim=c(0,1))

coverage_subsets_4_v4_r <- poisson_simulation_subsets4(100,250,4,0.05,R=40)
plot(coverage_subsets_4_v4_r, ylim=c(0,1))
