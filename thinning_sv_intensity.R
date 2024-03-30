library(spatstat)
library(latex2exp)
library(unikn)

# using thinning to obtain confidence intervals
# this time instead of basic bootstrap, use the sample variance scaled by (nr/d)

thinning_sv_intensity <- function(data, thinning_param, alpha, R=99, W=1, r=seq(0.0, 0.1), df) {
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
    # subprocess <- as.ppp(subprocess_df, W=square(W))
    i <- nrow(subprocess_df)/(W^2 * thinning_param)
    i_vals <- c(i_vals, i)
  }
  i_est <- npoints(data)/W^2
  # scalar <- npoints*thinning_param/(1-thinning_param)
  # scalar <- npoints*nrow(subprocess_df)/(npoints-nrow(subprocess_df))
  # scalar <- thinning_param
  scalar <- thinning_param/(1-thinning_param)
  # t <- qt((1-alpha/2), df = 1/thinning_param)
  t <- qt((1-alpha/2), df = R-1)
  lower_approx <- c()
  upper_approx <- c()
  
  lower_approx <- c(lower_approx, i_est - t*sqrt(var(i_vals)*scalar))
  upper_approx <- c(upper_approx, i_est + t*sqrt(var(i_vals)*scalar))

  return(cbind(lower_approx, upper_approx))
}



poisson_expanding_window_sv_intensity <- function(nsim, thinning_param, alpha, intensity, R=99, df=98) {
  
  # specifying window sizes over which to simulate
  Ws <- seq(0.3,3,0.2)
  cover <- rep(c(0),each=length(Ws))
  coverage <- cbind(Ws, cover)
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity, win=square(Ws[i]))
      confidences <- thinning_sv_intensity(p, thinning_param, alpha, W=Ws[i],r=c(0.0,0.1), R=R, df=df)
      if (confidences[1] <= intensity & intensity <= confidences[2]) {
        coverage[i,2] <- coverage[i,2] + 1/nsim
      }
    }
    print(coverage[i,])
  }
  # print(coverage)
  return(coverage)
}

intensity_cover_25 <- poisson_expanding_window_sv_intensity(2000, 0.25, 0.05, 250, 500)
intensity_cover_5 <- poisson_expanding_window_sv_intensity(2000, 0.5, 0.05, 250, 500)
intensity_cover_75 <- poisson_expanding_window_sv_intensity(2000, 0.75, 0.05, 250, 500)
save(intensity_cover_25, file="intensity_scaled_hom_pois_25.RData")
save(intensity_cover_5, file="intensity_scaled_hom_pois_5.RData")
save(intensity_cover_75, file="intensity_scaled_hom_pois_75.RData")

intensity_cover_25_r_over_d <- poisson_expanding_window_sv_intensity(2000, 0.25, 0.05, 250, 500)
intensity_cover_5_r_over_d <- poisson_expanding_window_sv_intensity(2000, 0.5, 0.05, 250, 500)
intensity_cover_75_r_over_d <- poisson_expanding_window_sv_intensity(2000, 0.75, 0.05, 250, 500)
save(intensity_cover_25_r_over_d, file="intensity_scaled_hom_pois_25_r_over_d.RData")
save(intensity_cover_5_r_over_d, file="intensity_scaled_hom_pois_5_r_over_d.RData")
save(intensity_cover_75_r_over_d, file="intensity_scaled_hom_pois_75_r_over_d.RData")

intensity_cover_25_est <- poisson_expanding_window_sv_intensity(2000, 0.25, 0.05, 50, 500)
intensity_cover_5_est <- poisson_expanding_window_sv_intensity(2000, 0.5, 0.05, 50, 500)
intensity_cover_75_est <- poisson_expanding_window_sv_intensity(2000, 0.75, 0.05, 50, 500)
save(intensity_cover_25_est, file="intensity_scaled_hom_pois_25_est.RData")
save(intensity_cover_5_est, file="intensity_scaled_hom_pois_5_est.RData")
save(intensity_cover_75_est, file="intensity_scaled_hom_pois_75_est.RData")

palette("default")   
par(mfrow=c(1,3))
plot(intensity_cover_25_est,type="l",lwd=1.5,ylim=c(0.5,1),xlab="",ylab="",font.main=2)
lines(intensity_bb_25_2,lwd=1.5, col=3)
title(mgp=c(2.5,0,0),xlab="Window Side Length",ylab="Confidence Interval Cover")
title(main=TeX(r'(p=0.25)'),line=1)
abline(h=0.95,col=2,lty=2,lwd=1.5)
plot(intensity_cover_5_est,type="l",lwd=1.5,ylim=c(0.5,1),xlab="",ylab="")
lines(intensity_bb_5_2,lwd=1.5, col=3)
title(mgp=c(2.5,0,0),xlab="Window Side Length",ylab="Confidence Interval Cover", )
title(main=TeX(r'(p=0.5)'),line=1)
abline(h=0.95,col=2,lty=2,lwd=1.5)
plot(intensity_cover_75_est,type="l",lwd=1.5,ylim=c(0.5,1),xlab="",ylab="")
lines(intensity_bb_75_2,lwd=1.5, col=3)
title(mgp=c(2.5,0,0),xlab="Window Side Length",ylab="Confidence Interval Cover", )
title(main=TeX(r'(p=0.75)'),line=1)
abline(h=0.95,col=2,lty=2,lwd=1.5)
mtext(font=2,text="Confidence Interval Covers For Estimating Intensity", side = 3, line = -2, outer = TRUE)
legend(1.8,0.55,c("Scaled CI","Naive CI"),col=c(1,3),lty=c(1,1),lwd=c(1.5,1.5))

