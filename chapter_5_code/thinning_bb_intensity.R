# scaling the thinned estimates before using basic bootstrap
# not working well

thinning_bb_intensity <- function(data, thinning_param, alpha, R=99, W=1, r=seq(0.0, 0.1), df) {
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
  i_vals <- c()
  for (i in 1:R) {
    unif <- runif(npoints, 0, 1)
    subprocess_df <- process_df[which(unif < thinning_param),]
    i <- nrow(subprocess_df)/(W^2*thinning_param)
    i_vals <- c(i_vals, i)
  }
  i_est <- npoints/W^2
  lower_approx <- c()
  upper_approx <- c()
  
  conf_interval <- unname(quantile(i_vals, probs=c(0.025,0.975), na.rm=TRUE))
  lower_approx <- c(lower_approx, 2*i_est - conf_interval[2])
  upper_approx <- c(upper_approx, 2*i_est - conf_interval[1])
  return(cbind(lower_approx, upper_approx))
}



poisson_expanding_window_bb_intensity <- function(nsim, thinning_param, alpha, intensity, R=99, df=98) {
  
  # specifying window sizes over which to simulate
  # Ws <- seq(0.5,5,0.5)
  Ws <- seq(0.5,2,0.25)
  cover <- rep(c(0),each=length(Ws))
  coverage <- cbind(Ws, cover)
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity, win=square(Ws[i]))
      confidences <- thinning_bb_intensity(p, thinning_param, alpha, W=Ws[i],r=c(0.0,0.1), R=R, df=df)
      if (confidences[1] <= intensity & intensity <= confidences[2]) {
        coverage[i,2] <- coverage[i,2] + 1/nsim
      }
    }
    print(coverage[i,])
  }
  # print(coverage)
  return(coverage)
}

intensity_bb_2_2 <- poisson_expanding_window_bb_intensity(1000,0.2,0.05,50,R=500)
intensity_bb_5_2 <- poisson_expanding_window_bb_intensity(1000,0.5,0.05,50,R=500)
intensity_bb_8_2 <- poisson_expanding_window_bb_intensity(1000,0.8,0.05,50,R=500)
save(intensity_bb_2_2, file="intensity_bb_2_2.RData")
save(intensity_bb_5_2, file="intensity_bb_5_2.RData")
save(intensity_bb_8_2, file="intensity_bb_8_2.RData")
