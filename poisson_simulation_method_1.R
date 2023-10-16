library(spatstat)

poisson_simulation_approximation <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    for (j in 1:length(r)) {
      data <- rpoispp(lambda)
      confidences <- approximation_method(data, nregions, alpha)
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_tiling <- function(nsim, lambda, nregions, alpha, R = 100) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    for (j in 1:length(r)) {
      data <- rpoispp(lambda)
      confidences <- tiling_method(data, nregions, alpha, R)
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_subsets <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    for (j in 1:length(r)) {
      data <- rpoispp(lambda)
      confidences <- subsets_method(data, nregions, alpha, R)
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_marked_point <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    for (j in 1:length(r)) {
      data <- rpoispp(lambda)
      confidences <- marked_point_method(data, nregions, alpha, R)
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

#regions4 <- poisson_simulation(1000, 250, 4, 0.05)
#regions16 <- poisson_simulation(1000, 250, 16, 0.05)
#regions64 <- poisson_simulation(1000, 250, 64, 0.05)

#plot(regions4["r"], regions4["cover"])
#lines(regions16$r, regions16$cover)
#lines(regions64$r, regions64$cover)