poisson_simulation_thinning_sv_p_squared <- function(nsim, lambda, thinning_param, alpha, R = 499, inhomogenous = FALSE, intensity_est = FALSE) {
  r <- seq(0.0, 0.14, 0.01)
  if (inhomogenous) {
    K_actual <- K_actual_inhom()
  } else {
    K_actual <- rep(c(pi),each=15) * r * r
  }
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  
  # scaler <- scaling_constant_inhom(thinning_param,intensity=intensity2, mean=TRUE)
  scaler <- 1/thinning_param^2
  for (i in 1:nsim) {
    tryCatch({print(paste0("Current simulation:",i))
      if (inhomogenous) {
        data <- rpoispp(intensity2)
      }else {
        data <- rpoispp(lambda)
      }
      confidences <- thinning_sample_var(data, thinning_param, alpha, R = R, inhomogenous = inhomogenous, scaler = scaler, intensity_est = intensity_est)
      for (j in 1:length(r)) {
        if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
          coverage[j,2] = coverage[j,2] + 1/nsim
        }
      }}, warning = function(w) { print("A simulation didn't work! Maybe a block is empty?") })
  }
  return(coverage)
}

thinning_sv_08_2 <- poisson_simulation_thinning_sv_p_squared(1000, 250, 0.8, 0.05)
thinning_sv_05_2 <- poisson_simulation_thinning_sv_p_squared(1000, 250, 0.5, 0.05)
thinning_sv_02_2 <- poisson_simulation_thinning_sv_p_squared(1000, 250, 0.2, 0.05)

sc_simulation_thinning_sv_p_squared <- function(nsim, lambda, thinning_param, alpha, R=99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- soft_core_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  scaler <- 1/thinning_param^2
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating softcore process
    data <- simulate_soft_core()
    
    confidences <- thinning_sample_var(data, thinning_param, alpha, R=99, scaler = scaler)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

sc_thinning_sv_08_2 <- sc_simulation_thinning_sv_p_squared(100, 250, 0.8, 0.05)
sc_thinning_sv_05_2 <- sc_simulation_thinning_sv_p_squared(100, 250, 0.5, 0.05)
sc_thinning_sv_02_2 <- sc_simulation_thinning_sv_p_squared(100, 250, 0.2, 0.05)

ns_simulation_thinning_sv_p_squared <- function(nsim, lambda, thinning_param, alpha, R=99) {
  r <- seq(0.0, 0.14, 0.01)
  K_actual <- neyman_scott_actual_k()
  cover <- rep(c(0),each=15)
  coverage <- cbind(r, cover)
  scaler <- 1/thinning_param^2
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    # generating matérn cluster field
    data <- rMatClust(25,0.1,10, win=owin(c(0,1),c(0,1)))
    confidences <- thinning_sample_var(data, thinning_param, alpha, R=R, scaler = scaler)
    print(confidences)
    for (j in 1:length(r)) {
      if (confidences[j,1] <= K_actual[j] & K_actual[j] <= confidences[j,2]) {
        coverage[j,2] <- coverage[j,2] + 1/nsim
      }
    }
  }
  return(coverage)
}

ns_thinning_08_2 <- ns_simulation_thinning_sv_p_squared(100, 250, 0.8, 0.05)
ns_thinning_05_2 <- ns_simulation_thinning_sv_p_squared(100, 250, 0.5, 0.05)
ns_thinning_02_2 <- ns_simulation_thinning_sv_p_squared(100, 250, 0.2, 0.05)


par(mfrow=c(1,3))

par(mgp=c(2,1,0))
plot(thinning_sv_08_2[-1,], type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, h", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(thinning_sv_05_2[-1,], type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(thinning_sv_02_2[-1,], type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
legend(0.01,0.25,c("p=0.8","p=0.5","p=0.2"),col=c(7,15,22),lwd=c(1.5,1.5,1.5),cex=1.3)
title("Homogenous Poisson",line=1,cex.main=1.5)

par(mgp=c(2,1,0))
plot(ns_thinning_08_2[-1,], type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, h", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(ns_thinning_05_2[-1,], type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(ns_thinning_02_2[-1,], type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
title("Neyman-Scott",line=1,cex.main=1.5)

par(mgp=c(2,1,0))
plot(sc_thinning_sv_08_2[-1,], type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, h", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(sc_thinning_sv_05_2[-1,], type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(sc_thinning_sv_02_2[-1,], type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
title("Soft Core",line=1,cex.main=1.5)

title("Studentized Cover for Ripley's K Using Approximate Scalar",line=-1.5,cex.main=1.5,outer=TRUE)

