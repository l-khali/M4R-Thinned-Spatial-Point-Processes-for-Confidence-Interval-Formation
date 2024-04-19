library(spatstat)

# using thinning to obtain confidence intervals
# this time instead of basic bootstrap, use the sample variance scaled by (nr/d)

thinning_sv_L <- function(data, thinning_param, alpha, R=99, W=1, r=seq(0.0, 0.1), df) {
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
  l_vals <- data.frame(r)
  for (i in 1:R) {
    unif <- runif(npoints, 0, 1)
    subprocess_df <- process_df[which(unif < thinning_param),]
    subprocess <- rthin(data, thinning_param)
    l <- Lest(subprocess, r=r, correction=c("isotropic"))
    l_vals <- cbind(l_vals, as.data.frame(l)["iso"])
  }
  L_est <- as.data.frame(Lest(data, r=r, correction=c("isotropic")))["iso"]
  # print(K_est[2,])
  # scalar <- npoints*thinning_param/(1-thinning_param)
  # scalar <- npoints*nrow(subprocess_df)/(npoints-nrow(subprocess_df))
  scalar <- thinning_param/((1-thinning_param))
  # t <- qt((1-alpha/2), df = 1/thinning_param)
  t <- qt((1-alpha/2), df = R-1)
  lower_approx <- c()
  upper_approx <- c()
  
  for (radius in 1:length(r)) {
    lower_approx <- c(lower_approx, L_est[radius,] - t*sqrt(var(as.numeric(l_vals[radius,-1]))*scalar))
    upper_approx <- c(upper_approx, L_est[radius,] + t*sqrt(var(as.numeric(l_vals[radius,-1]))*scalar))
  }
  return(cbind(lower_approx, upper_approx))
}



poisson_expanding_window_sv <- function(nsim, thinning_param, alpha, intensity, R=99, df=98) {
  
  # specifying window sizes over which to simulate
  Ws <- seq(0.5,2,0.25)
  cover <- rep(c(0),each=length(Ws))
  coverage <- cbind(Ws, cover)
  # indexing radii to select radius of 0.1
  j <- 2
  # actual Ripley's K for Thomas process
  L_actual <- 0.1
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity, win=square(Ws[i]))
      confidences <- thinning_sv_L(p, thinning_param, alpha, W=Ws[i],r=c(0.0,0.1), R=R, df=df)
      if (!is.na(confidences[j,1]) & !is.na(confidences[j,1])) {
        if (confidences[j,1] <= L_actual & L_actual <= confidences[j,2]) {
          coverage[i,2] <- coverage[i,2] + 1/nsim
        }
      }
    }
    print(coverage[i,])
  }
  # print(coverage)
  return(coverage)
}

