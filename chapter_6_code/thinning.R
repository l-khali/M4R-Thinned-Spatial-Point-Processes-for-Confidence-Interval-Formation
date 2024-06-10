library(spatstat)
library(pracma)
library(latex2exp)
library(rstudioapi)
library(docstring)
#' 
#' # thinning <- function(data, thinning_param, alpha, R=99) {
#'   #' Implement thinning method to obtain (1-alpha)*100% confidence intervals for
#'   #' a range of radii. The radii used are a sequence from 0.01 to 0.14 with a
#'   #' step of 0.01 as used in Loh and Stein (2004). Uses the basic bootstrap to
#'   #' calculate the confidence intervals.
#'   #' 
#'   #' thinning_param: determines what proportion of the samples are retained
#'   #' during thinning
#'   #' alpha: confidence level
#'   #' R: the number of thinned samples to take in order to form confidence 
#'   #' intervals
#'   #' 
#' 
#' # thinning <- function(data, thinning_param, alpha, R=99, W=1) {
#' #   # confidence intervals calculated using quantiles of samples
#' #   process_df <- as.data.frame(data)
#' #   npoints <- nrow(process_df)
#' #   r <- seq(0.0, 0.14, 0.01)
#' #   k_vals <- data.frame(r)
#' #   for (i in 1:R) {
#' #     unif <- runif(npoints, 0, 1)
#' #     subprocess_df <- process_df[which(unif < thinning_param),]
#' #     subprocess <- as.ppp(subprocess_df, owin(c(0,W),c(0,W)))
#' #     k <- Kest(subprocess, r=r, correction=c("isotropic"))
#' #     k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
#' #   }
#' #   K_est <- as.data.frame(Kest(data, r=r, correction=c("isotropic")))["iso"]
#' #   lower_approx <- c()
#' #   upper_approx <- c()
#' #   
#' #   for (radius in 1:15) {
#' #     sorted_k_vals <- sort(as.numeric(k_vals[radius,-1]))
#' #     lower <- sorted_k_vals[(R+1)*(1-alpha/2)]
#' #     upper <- sorted_k_vals[(R+1)*(alpha/2)]
#' #     lower_approx <- c(lower_approx, 2*K_est[radius,] - lower)
#' #     upper_approx <- c(upper_approx, 2*K_est[radius,] - upper)
#' #   }
#' #   return(cbind(lower_approx, upper_approx)) 
#' # }
#' 
#' thinning <- function(data, thinning_param, alpha, R=99, W=1, r=seq(0.0, 0.14, 0.01)) {
#'   #' Implement thinning method to obtain (1-alpha)*100% confidence intervals for
#'   #' a range of radii. The radii used are a sequence from 0.01 to 0.14 with a
#'   #' step of 0.01 as used in Loh and Stein (2004). Uses the basic bootstrap to
#'   #' calculate the confidence intervals.
#'   #' 
#'   #' thinning_param: determines what proportion of the samples are retained
#'   #' during thinning
#'   #' alpha: confidence level
#'   #' R: the number of thinned samples to take in order to form confidence 
#'   #' intervals
#'   #' 
#'   
#'   # confidence intervals calculated using quantiles of samples
#'   process_df <- as.data.frame(data)
#'   npoints <- nrow(process_df)
#'   # r <- seq(0.0, 0.14, 0.01)
#'   k_vals <- data.frame(r)
#'   for (i in 1:R) {
#'     unif <- runif(npoints, 0, 1)
#'     subprocess_df <- process_df[which(unif < thinning_param),]
#'     subprocess <- as.ppp(subprocess_df, owin(c(0,W),c(0,W)))
#'     k <- Kest(subprocess, r=r, correction=c("isotropic"))
#'     k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
#'   }
#'   K_est <- as.data.frame(Kest(data, r=r, correction=c("isotropic")))["iso"]
#'   lower_approx <- c()
#'   upper_approx <- c()
#' 
#'   for (radius in 1:length(r)) {
#'     sorted_k_vals <- sort(as.numeric(k_vals[radius,-1]))
#'     lower <- sorted_k_vals[(R+1)*(1-alpha/2)]
#'     upper <- sorted_k_vals[(R+1)*(alpha/2)]
#'     lower_approx <- c(lower_approx, 2*K_est[radius,] - lower)
#'     upper_approx <- c(upper_approx, 2*K_est[radius,] - upper)
#'   }
#'   print(upper_approx)
#'   return(cbind(lower_approx, upper_approx))
#' }



thinning <- function(data, thinning_param, alpha, R=500, W=1, r=seq(0.0, 0.14, 0.01)) {
  #' Implement thinning method to obtain (1-alpha)*100% confidence intervals for
  #' a range of radii. The radii used are a sequence from 0.01 to 0.14 with a
  #' step of 0.01 as used in Loh and Stein (2004). Uses the basic bootstrap to
  #' calculate the confidence intervals.
  #' 
  #' thinning_param: determines what proportion of the samples are retained
  #' during thinning
  #' alpha: confidence level
  #' R: the number of thinned samples to take in order to form confidence 
  #' intervals
  #' 
  
  # confidence intervals calculated using quantiles of samples
  process_df <- as.data.frame(data)
  npoints <- nrow(process_df)
  # r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  for (i in 1:R) {
    unif <- runif(npoints, 0, 1)
    subprocess_df <- process_df[which(unif < thinning_param),]
    subprocess <- as.ppp(subprocess_df, owin(c(0,W),c(0,W)))
    k <- Kest(subprocess, r=r, correction=c("isotropic"))
    k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
  }
  K_est <- as.data.frame(Kest(data, r=r, correction=c("isotropic")))["iso"]
  print(length(k_vals))
  lower_approx <- c()
  upper_approx <- c()

  for (radius in 1:length(r)) {
    sorted_k_vals <- sort(as.numeric(k_vals[radius,-1]))
    print(sorted_k_vals)
    print(sorted_k_vals[floor((R+1)*(1-alpha/2))])
    print(sorted_k_vals[ceiling((R+1)*(1-alpha/2))])
    lower <- (sorted_k_vals[floor((R+1)*(1-alpha/2))] + sorted_k_vals[ceiling((R+1)*(1-alpha/2))])/2
    upper <- (sorted_k_vals[floor((R+1)*(alpha/2))] + sorted_k_vals[ceiling((R+1)*(alpha/2))]) / 2
    print(paste0("lower",lower))
    lower_approx <- c(lower_approx, 2*K_est[radius,] - lower)
    upper_approx <- c(upper_approx, 2*K_est[radius,] - upper)
  }
  return(cbind(lower_approx, upper_approx))
}




thinning_1 <- poisson_simulation_thinning(500, 250, 1, 0.05)
plot(thinning_1[-1,], type="l", lty=1, ylim=c(0,1), main="Poisson: thinning (1.0)")
thinning_09 <- poisson_simulation_thinning(1000, 250, 0.9, 0.05)
plot(thinning_09[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.9)")
thinning_08 <- poisson_simulation_thinning(1000, 250, 0.8, 0.05)
plot(thinning_08[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.8)")
thinning_07 <- poisson_simulation_thinning(500, 250, 0.7, 0.05)
plot(thinning_07[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.7)")
thinning_06 <- poisson_simulation_thinning(500, 250, 0.6, 0.05)
plot(thinning_06[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.6)")
thinning_05 <- poisson_simulation_thinning(1000, 250, 0.5, 0.05)
plot(thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.5)")
thinning_04 <- poisson_simulation_thinning(500, 250, 0.4, 0.05)
plot(thinning_04[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.4)")
thinning_03 <- poisson_simulation_thinning(500, 250, 0.3, 0.05)
plot(thinning_03[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.3)")
thinning_02 <- poisson_simulation_thinning(1000, 250, 0.2, 0.05)
plot(thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.2)")
thinning_01 <- poisson_simulation_thinning(500, 250, 0.1, 0.05)
plot(thinning_01[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.1)")

par(mgp=c(2,1,0))  
plot(thinning_09[-1,], type="l", lty=1, ylim=c(0.5,1), col=1, xlab="Radius", ylab="Confidence Interval Cover")
lines(thinning_08[-1,], type="l", lty=1, ylim=c(0.5,1), col=2)
lines(thinning_07[-1,], type="l", lty=1, ylim=c(0.5,1), col=3)
lines(thinning_06[-1,], type="l", lty=1, ylim=c(0.5,1), col=4)
lines(thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=5)
lines(thinning_04[-1,], type="l", lty=1, ylim=c(0.5,1), col=6)
lines(thinning_03[-1,], type="l", lty=1, ylim=c(0.5,1), col=7)
lines(thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=8)
lines(thinning_01[-1,], type="l", lty=1, ylim=c(0.5,1), col=9)
legend(0.11, 0.77, legend=c("p=0.9 ", "p=0.8", "p=0.7", "p=0.6", "p=0.5", "p=0.4", "p=0.3", "p=0.2", "p=0.1"),
       col=c(1,2,3,4,5,6,7,8,9), lty =1, cex=0.8)
abline(h=0.95,lty=2)
title("Cover on Homogenous Poisson Process", line = 0.3)

thinning_08_2 <- thinning_08[-1,]
thinning_05_2 <- thinning_05[-1,]
thinning_02_2 <- thinning_02[-1,]
save(thinning_08_2, file="thinning_08_naive")
save(thinning_05_2, file="thinning_05_naive")
save(thinning_02_2, file="thinning_02_naive")

thinning_sample_var <- function(data, thinning_param, alpha, R=99, inhomogenous=FALSE, scaler = 1, intensity_est = FALSE) {
  # confidence intervals calculated using sample variance
  process_df <- as.data.frame(data)
  npoints <- nrow(process_df)
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  for (i in 1:R) {
    unif <- runif(npoints, 0, 1)
    subprocess_df <- process_df[which(unif < thinning_param),]
    subprocess <- as.ppp(subprocess_df, owin(c(0,1),c(0,1)))
    if (inhomogenous) {
      if (intensity_est) {
        k <- Kinhom(subprocess, r=r, correction=c("isotropic"))
      } else {
        k <- Kinhom(subprocess, lambda=intensity2, r=r, correction=c("isotropic"))
      }
    } else {
      k <- Kest(subprocess, r=r, correction=c("isotropic"))
    }
    k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
  }
  
  # scaler <- scaling_constant_inhom(thinning_param, 250, mean=TRUE)
  if (inhomogenous) {
    K_est <- as.data.frame(Kinhom(data, lambda=intensity2, r=r, correction=c("isotropic")))["iso"]
  } else {
    K_est <- as.data.frame(Kest(data, r=r, correction=c("isotropic")))["iso"]
  }
  t <- qt((1-alpha/2), df = R - 1)
  
  lower_approx <- c()
  upper_approx <- c()
  
  for (radius in 1:15) {
    sample_variance <- var(as.numeric(k_vals[radius,-1]))
    lower_approx <- c(lower_approx, K_est[radius,] - t * sqrt(sample_variance/scaler))
    upper_approx <- c(upper_approx, K_est[radius,] + t * sqrt(sample_variance/scaler))
  }
  return(cbind(lower_approx, upper_approx)) 
}



