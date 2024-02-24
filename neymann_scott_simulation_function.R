library(spatstat)
source("method_2_bootstrapping_spatial_point_patterns.R")
source("method_4_bootstrapping_spatial_processes.R")

# uniformly generated samples on a disc
#nclust <-  function(x0, y0, radius) {
#  n <- rpois(1,10)
#  return(runifdisc(n, radius, centre=c(x0, y0)))
#}

neyman_scott_actual_k <- function(r=seq(0.0, 0.14, 0.01), lambda=25, R=0.1) {
  K_actual <- c()
  for (radius in 1:length(r)) {
    rad <- r[radius]
    if (rad > 2*R) {
      K_actual <- c(K_actual, pi*rad^2 + 1/lambda)
    } else {
      z <- rad/(2*R)
      k <- pi*(rad^2) + (1/lambda)*(2 + (1/pi)*( ((8*(z^2) - 4)*acos(z)) - 2*asin(z) + 4*z*sqrt((1 - z^2)^3) - 6*z*sqrt(1 - (z^2))))
      K_actual <- c(K_actual, k)
    }
  }
  return(K_actual)
}

ns_simulation_approximation <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- neyman_scott_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
  
    # generating matérn cluster field
    data <- rMatClust(25,0.1,10, win=owin(c(0,1),c(0,1)))
    
    confidences <- approximation_method_3(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}


ns_simulation_tiling <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- neyman_scott_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    # generating matérn cluster field
    data <- rMatClust(25,0.1,10, win=owin(c(0,1),c(0,1)))
    confidences <- tiling_method(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

ns_simulation_subsets5 <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- neyman_scott_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    # generating matérn cluster field
    data <- rMatClust(25,0.1,10, win=owin(c(0,1),c(0,1)))
    confidences <- subsets_method5(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

ns_simulation_marked <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- neyman_scott_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    # generating matérn cluster field
    data <- rMatClust(25,0.1,10, win=owin(c(0,1),c(0,1)))
    confidences <- marked_point_method(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

ns_simulation_marked2 <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- neyman_scott_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    # generating matérn cluster field
    data <- rMatClust(25,0.1,10, win=owin(c(0,1),c(0,1)))
    confidences <- marked_point_method2(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

ns_simulation_thinning <- function(nsim, lambda, thinning_param, alpha, R=99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- neyman_scott_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    # generating matérn cluster field
    data <- rMatClust(25,0.1,10, win=owin(c(0,1),c(0,1)))
    confidences <- thinning(data, thinning_param, alpha, R=R)
    print(confidences)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

ns_simulation_thinning_sv <- function(nsim, lambda, thinning_param, alpha, R=99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- neyman_scott_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    # generating matérn cluster field
    data <- rMatClust(25,0.1,10, win=owin(c(0,1),c(0,1)))
    confidences <- thinning_sample_var(data, thinning_param, alpha, R=R)
    print(confidences)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}
