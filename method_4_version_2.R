source("poisson_simulation_method_1.R")
library(spatstat)

marked_point_method2 <- function(data, N, alpha, R=99) {
  npoints <- nrow(as.data.frame(data))
  r <- t(as.data.frame(replicate(npoints,seq(0.0, 0.14, 0.01))))
  weights <- edge.Ripley(data,pairdist(data))
  a <- 1/N
  increment <- 1/sqrt(N)
  window <- owin(c(0, increment), c(0,increment))
  r_vector <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r_vector)
  
  marks <- as.data.frame(matrix(0, ncol = 15, nrow = npoints))
  datadf <- as.data.frame(data)
  
  for (point in 1:(nrow(datadf)-1)) {
    x1 <- datadf[point,1]
    y1 <- datadf[point,2]
    for (other_point in (point+1):(nrow(datadf))) {
      x2 <- datadf[other_point,1]
      y2 <- datadf[other_point,2]
      for (radius_index in 1:15) {
        current_radius <- seq(0.0, 0.14, 0.01)[radius_index]
        if ((x1-x2)**2+(y1-y2)**2 <= current_radius**2) {
          marks[point,radius_index:15] <- marks[point,radius_index:15] + rep(weights[point, other_point], each=length(radius_index:15))
          marks[other_point,radius_index:15] <- marks[other_point,radius_index:15] + rep(weights[other_point, point], each=length(radius_index:15))
          break
        }
      }
    }
  }
  marks(data) <- marks
  
  # assuming that confence intervals formed in the same way as subsets method
  for (sim in 1:R) {
    n_sum <- 0
    k_vals_current <- data.frame(rep(c(0), each=15))
    for (i in 1:N) {
      xstart <- runif(1,0,1)
      ystart <- runif(1,0,1)
      xend <- xstart + 1/sqrt(N)
      yend <- ystart + 1/sqrt(N)
      
      # using modulus so that points that "wrap around" are next to each other
      if (xend <= 1 & yend <= 1) {
        subregion <- as.data.frame(subset(data, (xstart <= x & x < xend) & (ystart <= y & y < yend)))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
      } else if (xend > 1 & yend <= 1) {
        subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y & y < yend)))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart)))
      } else if (xend <= 1 & yend > 1) {
        subregion <- as.data.frame(subset(data, (xstart <= x & x < xend) & (ystart <= y | y < yend %% 1)))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart)))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
      } else {
        subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y | y < yend %% 1)))
        subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1))
        subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1))
      }
      
      n <- nrow(subregion)
      n_sum <- n_sum + n
      subregion <- as.ppp(subregion, window)
      #print(paste0("marks:",marks(subregion)))
      k <- colSums(marks(subregion))
      #print(paste0("k:",k))
      k_vals_current <- cbind(k_vals_current, as.numeric(k))
    }
    k_vals <- cbind(k_vals, rowSums(k_vals_current) / (n_sum * (n_sum - 1)))
  }
  
  K_est <- as.data.frame(Kest(data, r = r_vector, correction=c("isotropic")))["iso"]
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

#marked4 <- poisson_simulation_marked_point2(1000, 250, 4, 0.05)
#plot(marked4[-1,], ylim=c(0,1), type="l", lty=3, main="Poisson: marked point")
#marked16 <- poisson_simulation_marked_point2(1000, 250, 16, 0.05)
#lines(marked16[-1,], ylim=c(0,1), type="l", lty=2, main="Poisson: marked point")
#marked64 <- poisson_simulation_marked_point2(1000, 250, 64, 0.05)
#lines(marked64[-1,], ylim=c(0,1), type="l", lty=1, main="Poisson: marked point")

