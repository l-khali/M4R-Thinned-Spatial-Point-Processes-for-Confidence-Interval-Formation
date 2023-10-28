library(spatstat)

soft_core_pdf <- function(r) {
  if (0<=r & r<=0.05) {
    return(800*r)
  } else {
    return(0)
  }
}

simulate_soft_core <- function() {
  # Generate a Poisson point process
  poisson_pp <- rpoispp(500, win=owin(c(-0.1,1.1),c(-0.1,1.1)))
  npoints <- nrow(as.data.frame(poisson_pp))
  
  # Assign random markers to each point
  poisson_pp$marks <- runif(npoints, 0, 1)
  
  # Create a distance matrix between points
  dist_matrix <- crossdist(poisson_pp, poisson_pp)
  
  # Identify points with the minimum mark within the interaction radius
  i <- 1
  while (i <= npoints(poisson_pp)) {
    # if we have iterated through all points (because we are deleting as we go along), stop
    if (nrow(dist_matrix) < i) {
      break
    }
    # using inverse sampling method:
    neighbors <- which(dist_matrix[i,] < sqrt(runif(1,0,0.05)/20))
    if (length(neighbors) > 0) {
      #min_neighbor <- neighbors[which.min(poisson_pp$marks[neighbors])]
      #poisson_pp <- poisson_pp[-min_neighbor, ]
      #dist_matrix <- dist_matrix[-min_neighbor,]
      #dist_matrix <- dist_matrix[,-min_neighbor]
      if (min(poisson_pp$marks[neighbors]) < poisson_pp$marks[i]) {
        poisson_pp <- poisson_pp[-i,]
        dist_matrix <- dist_matrix[-i,]
        dist_matrix <- dist_matrix[,-i]
      } else {
      i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  
  # removing marks
  poisson_pp$marks <- rep(c(0), each=nrow(as.data.frame(poisson_pp)))
  
  # selecting subset in the required unit square
  soft_core <- subset(poisson_pp, 0 <= x & x <= 1 & 0 <= y & y <= 1)
  Window(soft_core) <- owin(c(0,1),c(0,1))
  
  return(soft_core)
}


soft_core_actual_k <- function() {
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  for (i in 1:10) {
    sim <- simulate_soft_core()
    k_vals <- cbind(k_vals, Kest(sim, r = r, correction=c("isotropic"))$iso)
  }
  k_actual <- c()
  for (radius in 1:length(r)) {
    k_actual <- c(k_actual, mean(as.numeric(k_vals[radius,-1])))
  }
  return(k_actual)
}

sc_simulation_approximation <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- soft_core_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating softcore process
    data <- simulate_soft_core()
    
    confidences <- approximation_method_3(data, nregions, alpha)
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
  K_actual <- soft_core_actual_k()
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

sc_simulation_subsets <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- soft_core_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating softcore process
    data <- simulate_soft_core()
    
    confidences <- subsets_method5(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

sc_simulation_marked <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- soft_core_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating softcore process
    data <- simulate_soft_core()
    
    confidences <- marked_point_method2(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

sc_simulation_thinning <- function(nsim, lambda, thinning_param, alpha, R=99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- soft_core_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating softcore process
    data <- simulate_soft_core()
    
    confidences <- thinning(data, thinning_param, alpha, R=R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

sc_simulation_thinning_sv <- function(nsim, lambda, thinning_param, alpha, R=99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- soft_core_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating softcore process
    data <- simulate_soft_core()
    
    confidences <- thinning_sample_var(data, thinning_param, alpha, R=99)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}