# scaling the thinned estimates before using basic bootstrap
# not working well

thinning_bb <- function(data, thinning_param, alpha, R=99, W=1, r=seq(0.0, 0.1), df) {
  #' Implement thinning method to obtain (1-alpha)*100% confidence intervals for
  #' a range of radii. The radii used are a sequence from 0.01 to 0.14 with a
  #' step of 0.01 as used in Loh and Stein (2004). Uses the smaple variance of 
  #' estimates to calculate the confidence intervals.
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
    subprocess <- as.ppp(subprocess_df, W=square(W))
    k <- Kest(subprocess, r=r, correction=c("isotropic"))
    k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
  }
  K_est <- as.data.frame(Kest(data, r=r, correction=c("isotropic")))["iso"]
  # scalar <- npoints*thinning_param/(1-thinning_param)
  # scalar <- npoints*nrow(subprocess_df)/(npoints-nrow(subprocess_df))
  scalar <- thinning_param
  lower_approx <- c()
  upper_approx <- c()
  
  for (radius in 1:length(r)) {
    scaled_k_vals <- k_vals[radius,-1] * scalar
    conf_interval <- unname(quantile(as.numeric(scaled_k_vals), probs=c(0.025,0.975), na.rm=TRUE))
    lower_approx <- c(lower_approx, 2*K_est[radius,] - conf_interval[2])
    upper_approx <- c(upper_approx, 2*K_est[radius,] - conf_interval[1])
  }
  return(cbind(lower_approx, upper_approx))
}



poisson_expanding_window_bb <- function(nsim, thinning_param, alpha, intensity, R=99, df=98) {
  
  # specifying window sizes over which to simulate
  Ws <- seq(0.5,5,0.2)
  cover <- rep(c(0),each=length(Ws))
  coverage <- cbind(Ws, cover)
  # indexing radii to select radius of 0.1
  j <- 2
  # actual Ripley's K for Thomas process
  K_actual <- pi*0.1*0.1
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity, win=square(Ws[i]))
      confidences <- thinning_bb(p, thinning_param, alpha, W=Ws[i],r=c(0.0,0.1), R=R, df=df)
      print(confidences)
      if (confidences[j,1] <= K_actual & K_actual <= confidences[j,2]) {
        coverage[i,2] <- coverage[i,2] + 1/nsim
      }
    }
    print(coverage[i,])
  }
  # print(coverage)
  return(coverage)
}

cover_scaled_bb_9 <- poisson_expanding_window_bb(10,0.9,0.05,250,R=50)
plot(cover_scaled_bb_9, main="Scaled Conf Interval, df=R-1", type="l")

cover_scaled_bb_5 <- poisson_expanding_window_bb(10,0.5,0.05,250,R=50)
plot(cover_scaled_bb_5, main="Scaled Conf Interval, df=R-1", type="l")

cover_scaled_bb_1 <- poisson_expanding_window_bb(10,0.1,0.05,250,R=50)
plot(cover_scaled_bb_1, main="Scaled Conf Interval, df=R-1", type="l")
