library(spatstat)

intensity2 <- function(x, y) {
  return(500*x)
}

intensity <- function(x,y) {
  return(240 + 20*x)
}

#poisson_simulation_approximation <- function(nsim, lambda, nregions, alpha) {
#   r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- approximation_method(data, nregions, alpha)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }

#poisson_simulation_approximation2 <- function(nsim, lambda, nregions, alpha) {
#   r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- approximation_method_2(data, nregions, alpha)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }

poisson_simulation_approximation3 <- function(nsim, lambda, nregions, alpha, inhomogenous = FALSE) {
  r <- seq(0.0, 0.14, 0.01)
  if (inhomogenous) {
    K_actual <- K_actual_inhom()
  } else {
    K_actual <- rep(c(pi),each=15) * r * r
  }
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    if (inhomogenous) {
      data <- rpoispp(intensity2)
    }else {
      data <- rpoispp(lambda)
    }
    confidences <- approximation_method_3(data, nregions, alpha)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

poisson_simulation_tiling <- function(nsim, lambda, nregions, alpha, R = 100, inhomogenous = FALSE) {
  r <- seq(0.0, 0.14, 0.01)
  if (inhomogenous) {
    K_actual <- K_actual_inhom()
  } else {
    K_actual <- rep(c(pi),each=15) * r * r
  }
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    if (inhomogenous) {
      data <- rpoispp(intensity2)
    }else {
      data <- rpoispp(lambda)
    }
    confidences <- tiling_method(data, nregions, alpha, R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}
# 
# poisson_simulation_subsets <- function(nsim, lambda, nregions, alpha, R = 99) {
#   r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- subsets_method(data, nregions, alpha, R)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }
# 
# poisson_simulation_subsets2 <- function(nsim, lambda, nregions, alpha, R = 99) {
#   r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- subsets_method2(data, nregions, alpha, R)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }
# 
# poisson_simulation_subsets3 <- function(nsim, lambda, nregions, alpha, R = 99) {
#   r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- subsets_method3(data, nregions, alpha, R)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }
# 
# poisson_simulation_subsets4 <- function(nsim, lambda, nregions, alpha, R = 99) {
#   r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- subsets_method4(data, nregions, alpha, R=R)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }

poisson_simulation_subsets5 <- function(nsim, lambda, nregions, alpha, R = 99, inhomogenous = FALSE) {
  r <- seq(0.0, 0.14, 0.01)
  if (inhomogenous) {
    K_actual <- K_actual_inhom()
  } else {
    K_actual <- rep(c(pi),each=15) * r * r
  }
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    if (inhomogenous) {
      data <- rpoispp(intensity2)
    }else {
      data <- rpoispp(lambda)
    }
    confidences <- subsets_method5(data, nregions, alpha, R=R)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] = coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}
# 
# poisson_simulation_subsets6 <- function(nsim, lambda, nregions, alpha, R = 99) {
#   r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- subsets_method6(data, nregions, alpha, R=R)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }
# 
# poisson_simulation_subsets7 <- function(nsim, lambda, nregions, alpha, R = 99) {
#   r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- subsets_method7(data, nregions, alpha, R=R)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }

# poisson_simulation_marked_point <- function(nsim, lambda, nregions, alpha, R = 99) {
  # r <- seq(0.0, 0.14, 0.01)
#   K_actual <- rep(c(pi),each=15) * r * r
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     data <- rpoispp(lambda)
#     confidences <- marked_point_method(data, nregions, alpha, R)
#     for (j in 1:length(r)) {
#       if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }

poisson_simulation_marked_point2 <- function(nsim, lambda, nregions, alpha, R = 99, inhomogenous=FALSE) {
  r <- seq(0.0, 0.14, 0.01)
  if (inhomogenous) {
    K_actual <- K_actual_inhom()
  } else {
    K_actual <- rep(c(pi),each=15) * r * r
  }
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    if (inhomogenous) {
      data <- rpoispp(intensity2)
    }else {
      data <- rpoispp(lambda)
    }
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

K_actual_inhom <- function() {
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  for (i in 1:1000) {
    sim <- rpoispp(intensity2)
    k_vals <- cbind(k_vals, Kinhom(sim, lambda=intensity2, r = r, correction=c("isotropic"))$iso)
  }
  k_actual <- c()
  for (radius in 1:length(r)) {
    k_actual <- c(k_actual, mean(as.numeric(k_vals[radius,-1])))
  }
  return(k_actual)
}

poisson_simulation_thinning_sv <- function(nsim, lambda, thinning_param, alpha, R = 99, inhomogenous = FALSE) {
  r <- seq(0.0, 0.14, 0.01)
  if (inhomogenous) {
    K_actual <- K_actual_inhom()
  } else {
    K_actual <- rep(c(pi),each=15) * r * r
  }
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  
  scaler <- scaling_constant_inhom(thinning_param,intensity=intensity2, mean=TRUE)
  
  for (i in 1:nsim) {
        tryCatch({print(paste0("Current simulation:",i))
        if (inhomogenous) {
          data <- rpoispp(intensity2)
        }else {
          data <- rpoispp(lambda)
        }
        confidences <- thinning_sample_var(data, thinning_param, alpha, R = R, inhomogenous = inhomogenous, scaler = scaler)
        for (j in 1:length(r)) {
          if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
            coverage[j,2] = coverage[j,2] + 1/nsim
          }
        }}, warning = function(w) { print("A simulation didn't work! Maybe a block is empty?") })
  }
  return(coverage)
}

# poisson_simulation_thinning_L <- function(nsim, lambda, thinning_param, alpha, R = 99, inhomogenous = FALSE) {
#   r <- seq(0.0, 0.14, 0.01)
#   if (inhomogenous) {
#     L_actual <- sqrt(K_actual_inhom()/pi)
#   } else {
#     L_actual <- sqrt((rep(c(pi),each=15) * r * r) / pi)
#   }
#   cover <- rep(c(0),each=15)
#   coverage <- cbind(r, cover)
#   
#   scaler <- scaling_constant_inhom(thinning_param,intensity=intensity2, mean=TRUE)
#   
#   for (i in 1:nsim) {
#     print(paste0("Current simulation:",i))
#     if (inhomogenous) {
#       data <- rpoispp(intensity2)
#     }else {
#       data <- rpoispp(lambda)
#     }
#     confidences <- thinning_sample_var(data, thinning_param, alpha, R = R, inhomogenous = inhomogenous, scaler = scaler)
#     confidences_L <- sqrt(confidences/pi)
#     for (j in 1:length(r)) {
#       if (confidences_L[j,1] <= L_actual[j] & L_actual[j] <= confidences_L[j,2]) {
#         coverage[j,2] = coverage[j,2] + 1/nsim
#       }
#     }
#   }
#   return(coverage)
# }
