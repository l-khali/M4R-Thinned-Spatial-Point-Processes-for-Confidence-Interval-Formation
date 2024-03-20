library(spatstat)

thinning_sv_thomas <- function(data, thinning_param, alpha, R=99, W=1, r=c(0.0, 0.05), df) {
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
  # print(K_est[2,])
  # scalar <- npoints*thinning_param/(1-thinning_param)
  # scalar <- npoints*nrow(subprocess_df)/(npoints-nrow(subprocess_df))
  # scalar <- thinning_param^2
  scalar <- scalars[thinning_param*10]
  # t <- qt((1-alpha/2), df = 1/thinning_param)
  t <- qt((1-alpha/2), df = R-1)
  lower_approx <- c()
  upper_approx <- c()
  
  for (radius in 1:length(r)) {
    lower_approx <- c(lower_approx, K_est[radius,] - t*sqrt(var(as.numeric(k_vals[radius,-1]))*scalar))
    upper_approx <- c(upper_approx, K_est[radius,] + t*sqrt(var(as.numeric(k_vals[radius,-1]))*scalar))
  }
  return(cbind(lower_approx, upper_approx))
}



poisson_expanding_window_sv_thomas <- function(nsim, thinning_param, alpha, parent_intensity, radius_of_clusters, daughter_intensity, R=99, df=98) {
  
  # specifying window sizes over which to simulate
  Ws <- seq(0.5,3,0.5)
  cover <- rep(c(0),each=length(Ws))
  coverage <- cbind(Ws, cover)
  # indexing radii to select radius of 0.1
  j <- 2
  # actual Ripley's K for Thomas process
  # K_actual <- neyman_scott_actual_k(r=c(0.0,0.05),parent_intensity, radius_of_clusters, daughter_intensity,1)
  K_actual <- pi*r[j]*r[j] + (1/parent_intensity)*(1-exp((-r[j]^2)/(4*radius_of_clusters^2)))
  
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      # p <- rpoispp(intensity, win=square(Ws[i]))
      p <- rThomas(parent_intensity, radius_of_clusters, daughter_intensity, win=owin(c(0,Ws[i]),c(0,Ws[i])))
      confidences <- thinning_sv_thomas(p, thinning_param, alpha, W=Ws[i],r=c(0.0,0.05), R=R, df=df)
      if (!is.na(confidences[j,1]) & !is.na(confidences[j,1])) {
        if (confidences[j,1] <= K_actual & K_actual <= confidences[j,2]) {
          coverage[i,2] <- coverage[i,2] + 1/nsim
        }
      }
    }
    print(coverage[i,])
  }
  # print(coverage)
  return(coverage)
}

scaled_thomas_75 <- poisson_expanding_window_sv_thomas(500,0.75,0.05,25,0.1,10,R=100)
scaled_thomas_5 <- poisson_expanding_window_sv_thomas(500,0.5,0.05,25,0.1,10,R=100)
scaled_thomas_25 <- poisson_expanding_window_sv_thomas(500,0.25,0.05,25,0.1,10,R=100)

save(scaled_thomas_25,file="scaled_thomas_cover_25")
save(scaled_thomas_5,file="scaled_thomas_cover_5")
save(scaled_thomas_75,file="scaled_thomas_cover_75")

plot(scaled_thomas_25,type="l",col=3,ylim=c(0.5,1),lwd=1.5,main="Cover of Thomas Process With Scaling")
lines(scaled_thomas_5,col=4,lwd=1.5)
lines(scaled_thomas_75,col=7,lwd=1.5)
abline(h=0.95,col=2,lty=2)
legend(2.5,0.6,c("0.25","0.5","0.75"),c(3,4,7))


# par(mfrow=c(1,3))
# plot(cover25_2,type="l",lwd=0.5,ylim=c(0.5,1),xlab="",ylab="",font.main=2)
# lines(cover5_2,lwd=1.5,col=3)
# lines(cover5_2,lwd=1.5,col=3)


