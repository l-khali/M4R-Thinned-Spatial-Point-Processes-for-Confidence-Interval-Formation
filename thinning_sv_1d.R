library(spatstat)

thinning_sv_mean <- function(data, thinning_param, alpha, R=99, df=98) {
  #' Implement thinning method to obtain (1-alpha)*100% confidence intervals for
  #' the mean of a poisson distribution. Uses the smaple variance of 
  #' estimates to calculate the confidence intervals.
  #' 
  #' thinning_param: determines what proportion of the samples are retained
  #' during thinning
  #' alpha: confidence level
  #' R: the number of thinned samples to take in order to form confidence 
  #' intervals
  #' 
  
  # confidence intervals calculated using quantiles of samples
  npoints <- length(data)
  thinned_means <- c()
  for (i in 1:R) {
    unif <- runif(npoints, 0, 1)
    subprocess_df <- data[which(unif < thinning_param)]
    if (length(subprocess_df) < 1) {next}
    thinned_means <- cbind(thinned_means, mean(subprocess_df))
  }

  mean_est <- mean(data)
  # scalar <- npoints*length(subprocess_df)/(npoints-length(subprocess_df))
  # scalar <- length(subprocess_df)/(npoints-length(subprocess_df))
  scalar <- 1
  # t <- qt((1-alpha/2), df = npoints/(npoints-length(subprocess_df)))
  t <- qt((1-alpha/2), df = 1/(1-thinning_param))

  lower_approx <- mean(as.numeric(thinned_means)) - t*sqrt(var(as.numeric(thinned_means))/scalar)
  upper_approx <- mean(as.numeric(thinned_means)) + t*sqrt(var(as.numeric(thinned_means))/scalar)
  # print(t*sqrt(var(as.numeric(thinned_means))/scalar))

  return(c(lower_approx,upper_approx))
}



poisson_sv_1d <- function(nsim, thinning_param, alpha, intensity, R=99, df=98) {
  
  # specifying ns over which to simulate
  ns <- c(1000,2000,3000)
  cover <- rep(c(0),each=length(ns))
  coverage <- cbind(ns, cover)
  
  for (i in 1:length(ns)) {
    print(paste0("Current n:", ns[i]))
    for (sim in 1:nsim) {
      p <- rpois(n, intensity)
      confidences <- thinning_sv_mean(p, thinning_param, alpha, R=R, df=df)
      if (confidences[1] <= intensity & intensity <= confidences[2]) {
        coverage[i,2] <- coverage[i,2] + 1/nsim
      }
    }
    # print(coverage[i,])
  }
  # print(coverage)
  return(coverage)
}

cover <- poisson_sv_1d(100,0.2,0.05,10, R=500)
plot(cover)
