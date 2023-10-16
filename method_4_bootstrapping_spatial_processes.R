source("poisson_simulation_method_1.R")
library(spatstat)

marked_point_method <- function(data, N, alpha, R) {
  npoints <- nrow(as.data.frame(data))
  r <- t(as.data.frame(replicate(npoints,seq(0.0, 0.14, 0.01))))
  weights <- edge.Ripley(data,r)
  a <- 1/N
  increment <- 1/sqrt(N)
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
          marks[point,radius_index:15] <- marks[point,radius_index:15] + weights[other_point, radius_index:15]
          marks[other_point,radius_index:15] <- marks[other_point,radius_index:15] + weights[point, radius_index:15]
          break
        }
      }
    }
  }
  marks(data) <- marks
  
  # assuming that confence intervals formed in the same way as subsets method
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
    n <- nrow(as.data.frame(subregion))
    
    k <- colSums(marks(subregion)) * a / (n*(n-1))
    k_vals <- cbind(k_vals, as.numeric(k))
  }
  
  K_est <- as.data.frame(Kest(data, r = r_vector, correction=c("isotropic")))["iso"]
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

coverage4 <- poisson_simulation_marked_point(10, 250, 4, 0.05)
coverage16 <- poisson_simulation_marked_point(1000, 250, 16, 0.05)