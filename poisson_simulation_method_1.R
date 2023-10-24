library(spatstat)

poisson_simulation_approximation <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- approximation_method(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_approximation2 <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- approximation_method_2(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_approximation3 <- function(nsim, lambda, nregions, alpha) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- approximation_method_3(data, nregions, alpha)
    for (j in 1:length(r)) {
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
    data <- rpoispp(lambda)
    confidences <- tiling_method(data, nregions, alpha, R)
    for (j in 1:length(r)) {
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
    data <- rpoispp(lambda)
    confidences <- subsets_method(data, nregions, alpha, R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_subsets2 <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- subsets_method2(data, nregions, alpha, R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_subsets3 <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- subsets_method3(data, nregions, alpha, R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_subsets4 <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- subsets_method4(data, nregions, alpha, R=R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_subsets5 <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- subsets_method5(data, nregions, alpha, R=R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_subsets6 <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- subsets_method6(data, nregions, alpha, R=R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_subsets7 <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- subsets_method7(data, nregions, alpha, R=R)
    for (j in 1:length(r)) {
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
    data <- rpoispp(lambda)
    confidences <- marked_point_method(data, nregions, alpha, R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_marked_point2 <- function(nsim, lambda, nregions, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- marked_point_method2(data, nregions, alpha, R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_thinning <- function(nsim, lambda, thinning_param, alpha, R = 99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- rep(c(pi),each=15) * r * r
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    data <- rpoispp(lambda)
    confidences <- thinning(data, thinning_param, alpha, R = R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}