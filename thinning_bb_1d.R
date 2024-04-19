library(spatstat)

thinning_bb_mean <- function(data, thinning_param, alpha, R=99, df=98) {
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
    thinned_means <- cbind(thinned_means, var(subprocess_df))
  }
  
  mean_est <- var(data)
  sorted_means <- sort(thinned_means)
  lower <- (sorted_means[floor((R+1)*(1-alpha/2))] + 3*sorted_means[ceiling((R+1)*(1-alpha/2))])/4
  upper <- (sorted_means[floor((R+1)*(alpha/2))] + sorted_means[ceiling((R+1)*(alpha/2))]) / 2
  
  lower_approx <- 2*mean_est - lower
  upper_approx <- 2*mean_est - upper
  return(c(lower_approx,upper_approx))
}



poisson_bb_1d <- function(nsim, thinning_param, alpha, intensity, R=99, df=98) {
  
  # specifying ns over which to simulate
  ns <- seq(5,504,10)
  cover <- rep(c(0),each=length(ns))
  coverage <- cbind(ns, cover)
  
  for (i in 1:length(ns)) {
    print(paste0("Current n:", ns[i]))
    for (sim in 1:nsim) {
      p <- rpois(n, intensity)
      confidences <- thinning_bb_mean(p, thinning_param, alpha, R=R, df=df)
      if (confidences[1] <= intensity & intensity <= confidences[2]) {
        coverage[i,2] <- coverage[i,2] + 1/nsim
      }
    }
    # print(coverage[i,])
  }
  # print(coverage)
  return(coverage)
}

asymptotic_cover_8_bb <- poisson_bb_1d(1000,0.8,0.05,10,R=500)
plot(asymptotic_cover_8_bb, main="Asymptotic Cover Using Thinning, p=0.8", xlab="n", type="l")
abline(h=0.95,col=2,lty=2)

asymptotic_cover_5_bb <- poisson_bb_1d(1000,0.5,0.05,10,R=500)
plot(asymptotic_cover_5_bb, main="Asymptotic Cover Using Thinning, p=0.5", xlab="n", type="l")
abline(h=0.95,col=2,lty=2)

asymptotic_cover_2_bb <- poisson_bb_1d(1000,0.2,0.05,10,R=500)
plot(asymptotic_cover_2_bb, main="Asymptotic Cover Using Thinning, p=0.2", xlab="n", type="l")
abline(h=0.95,col=2,lty=2)

plot(asymptotic_cover_8_bb,ylim=c(0.6,1),type="l",col=7,lwd=1.5,xlab=TeX('Number of points, $n$'),ylab="Confidence Interval Cover")
lines(asymptotic_cover_5_bb,ylim=c(0.5,1),col=15,lwd=1.5)
lines(asymptotic_cover_2_bb,ylim=c(0.5,1),col=22,lwd=1.5)
abline(h=0.95, col=18,lwd=1.5,lty=2)
legend(750,0.85,c("p=0.8","p=0.5","p=0.2"),col=c(7,15,22),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),cex=0.9)
title("Basic Bootstrap CI Cover: Distribution",line=0.4)



asymptotic_cover_5_var <- poisson_bb_1d(100,0.5,0.01,10,R=500)
plot(asymptotic_cover_5_var)
abline(h=0.95,col=2)

asymptotic_cover_2_var <- poisson_bb_1d(100,0.2,0.05,10,R=500)
plot(asymptotic_cover_2_var)
abline(h=0.95,col=2)
