library(spatstat)
source("neymann_scott_simulation_function.R")
source("method_1_version_2.R")



simulate_soft_core <- function() {
  # Generate a Poisson point process
  poisson_pp <- rpoispp(500)
  npoints <- nrow(as.data.frame(poisson_pp))
  
  # Assign random markers to each point
  poisson_pp$marks <- runif(npoints, 0, 1)
  
  # Create a distance matrix between points
  dist_matrix <- crossdist(poisson_pp, poisson_pp)
  
  # Identify points with the minimum mark within the interaction radius
  for (i in 1:npoints) {
    # if we have iterated through all points, stop
    if (nrow(dist_matrix) < i) {
      break
    }
    neighbors <- which(dist_matrix[i,] < 800*runif(1,0,0.05))
    if (length(neighbors) > 1) {
      min_neighbor <- neighbors[which.min(poisson_pp$marks[neighbors])]
      poisson_pp <- poisson_pp[-min_neighbor, ]
      dist_matrix <- dist_matrix[-min_neighbor,]
      dist_matrix <- dist_matrix[,-min_neighbor]
    }
  }
  
  # removing marks
  poisson_pp$marks <- rep(c(0), each=nrow(as.data.frame(poisson_pp)))
  
  return(poisson_pp)
}

sc_simulation_approximation <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating softcore process
    data <- simulate_soft_core()
    
    confidences <- approximation_method_2(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}


sc_simulation_tiling <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating softcore process
    data <- simulate_soft_core()
    
    confidences <- tiling_method(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}