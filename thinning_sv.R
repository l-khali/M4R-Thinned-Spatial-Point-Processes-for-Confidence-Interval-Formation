library(spatstat)

# using thinning to obtain confidence intervals
# this time instead of basic bootstrap, use the sample variance scaled by (nr/d)

thinning_sv <- function(data, thinning_param, alpha, R=99, W=1, r=seq(0.0, 0.1), df) {
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
    subprocess <- rthin(data, thinning_param)
    retained <- npoints(subprocess)
    k <- Kest(subprocess, r=r, correction=c("isotropic"))
    k_vals_s <- as.data.frame(k)["iso"]
    # k_val <- k_vals_s*retained*(retained-1)/(thinning_param*npoints*(npoints-1))
    k_vals <- cbind(k_vals, k_vals_s)
  }
  K_est <- as.data.frame(Kest(data, r=r, correction=c("isotropic")))["iso"]
  # print(K_est[2,])
  # scalar <- npoints*thinning_param/(1-thinning_param)
  # scalar <- npoints*nrow(subprocess_df)/(npoints-nrow(subprocess_df))
  scalar <- thinning_param
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



poisson_expanding_window_sv <- function(nsim, thinning_param, alpha, intensity, R=99, df=98) {
  
  # specifying window sizes over which to simulate
  Ws <- seq(0.5,2,0.25)
  cover <- rep(c(0),each=length(Ws))
  coverage <- cbind(Ws, cover)
  # indexing radii to select radius of 0.1
  j <- 2
  # actual Ripley's K for Thomas process
  K_actual <- pi*0.1*0.1
  
  confidences_mat <- c()
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity, win=square(Ws[i]))
      confidences <- thinning_sv(p, thinning_param, alpha, W=Ws[i],r=c(0.0,0.1), R=R, df=df)
      if (i==length(Ws)) {
        confidences_mat <- rbind(confidences_mat, confidences[j,])
        print(confidences_mat)
      }
      if (!is.na(confidences[j,1]) & !is.na(confidences[j,1])) {
        if (confidences[j,1] <= K_actual & K_actual <= confidences[j,2]) {
          coverage[i,2] <- coverage[i,2] + 1/nsim
        }
      }
    }
    print(coverage[i,])
  }
  # print(coverage)
  return(list("coverage"=coverage,"CIs"=confidences_mat))
}

hom_pois_K_8 <- poisson_expanding_window_sv(1000,0.8,0.05,250,R=500)
save(hom_pois_K_8, file = "hom_pois_K_8.RData")

hom_pois_K_5 <- poisson_expanding_window_sv(1000,0.5,0.05,250,R=500)
save(hom_pois_K_5, file = "hom_pois_K_5.RData")

hom_pois_K_2 <- poisson_expanding_window_sv(1000,0.2,0.05,250,R=500)
save(hom_pois_K_2, file = "hom_pois_K_2.RData")

hom_pois_K_8_p_scaling <- poisson_expanding_window_sv(1000,0.8,0.05,250,R=500)
save(hom_pois_K_8_p_scaling, file = "hom_pois_K_8.RData")

hom_pois_K_5_scaling <- poisson_expanding_window_sv(1000,0.5,0.05,250,R=500)
save(hom_pois_K_5_scaling, file = "hom_pois_K_5.RData")

hom_pois_K_2_scaling <- poisson_expanding_window_sv(1000,0.2,0.05,250,R=500)
save(hom_pois_K_2_scaling, file = "hom_pois_K_2.RData")

plotCI(1:50, rep(pi*0.1*0.1,50),li=hom_pois_K_2$CIs[,1][1:50], ui=hom_pois_K_2$CIs[,2][1:50])


cover25 <- poisson_expanding_window_sv(2000,0.25,0.05,250,R=1000)
plot(cover25, type="l",ylim=c(0,1),xlab="",ylab="",lwd=2)
title(mgp=c(2.5,0,0), xlab="Window Side length",ylab="Cover", font.main= 4)
abline(h=0.95,col=2,lty=2,lwd=2)
save(cover25, file = "scaled_hom_pois_25.RData")

# cover3 <- poisson_expanding_window_sv(1000,0.3,0.05,250,R=500)
# plot(cover3, main="Scaled Conf Interval, df=R-1", type="l")

cover5 <- poisson_expanding_window_sv(100,0.5,0.05,250,R=100)
save(cover5, file = "scaled_hom_pois_5.RData")
# cover5 <- cover5_df2
plot(cover5, type="l",ylim=c(0,1),xlab="",ylab="",lwd=2)
title(mgp=c(2.5,0,0), xlab="Window Side length",ylab="Cover", font.main= 4)
abline(h=0.95,col=2,lty=2,lwd=2)

cover75 <- poisson_expanding_window_sv(2000,0.75,0.05,250,R=1000)
plot(cover75, type="l",ylim=c(0,1),xlab="",ylab="",lwd=2)
title(mgp=c(2.5,0,0), xlab="Window Side length",ylab="Cover", font.main= 4)
abline(h=0.95,col=2,lty=2,lwd=2)
save(cover75, file = "scaled_hom_pois_75.RData")

save(cover25,file="k_cover_scaled_25.RData")
save(cover5,file="k_cover_scaled_5.RData")
save(cover75,file="k_cover_scaled_75.RData")

cover25_2 <- poisson_expanding_window_sv(100,0.25,0.05,250,R=50)
cover5_2 <- poisson_expanding_window_sv(100,0.5,0.05,250,R=50)
cover75_2 <- poisson_expanding_window_sv(100,0.75,0.05,250,R=50)

plot(cover25_2,type="l",col=3,ylim=c(0.5,1),lwd=1.5,main="Cover of Ripleys K Hom Pois Process With Scaling")
lines(cover5_2,col=4,lwd=1.5)
lines(cover75_2,col=7,lwd=1.5)
abline(h=0.95,col=2,lty=2)
legend(2.5,0.6,c("0.25","0.5","0.75"),c(3,4,7))


# cover9 <- poisson_expanding_window_sv(1000,0.9,0.05,250,R=500)
# plot(cover9, main="Scaled Conf Interval, df=R-1", type="l")

# using df = N/d
# cover2_df2 <- poisson_expanding_window_sv(100,0.2,0.05,250,R=50)
# plot(cover2_df2, main="Scaled Conf Interval, df=N/d", type="l")
# 
# cover5_df2 <- poisson_expanding_window_sv(100,0.5,0.05,250,R=50) # actually with R-1 degrees of freedom
# plot(cover5_df2, main="Scaled Conf Interval, df=N/d", type="l")
# 
# cover9_df2 <- poisson_expanding_window_sv(100,0.9,0.05,250,R=50)
# plot(cover9_df2, main="Scaled Conf Interval, df=N/d", type="l")


# plotting final term which is not accounted for in the scaling
thinning_params <- seq(0.1,0.9,0.1)
diffs <- c()
for (thinning_param in thinning_params) {
  n <- 5 * 250
  np <- thinning_param * n
  used <- n*((2*area^2*beta)/(n^2))*(1+0.305*gamma+beta*(-1+0.0132*n*gamma))/thinning_param^2
  actual <- n*((2*area^2*beta)/(np^2))*(1+0.305*gamma+beta*(-1+0.0132*np*gamma))/thinning_param^2
  diffs <- c(diffs, abs(used - actual))
}

plot(diffs)






dfs <- seq(1,100)
covers <- cbind(seq(0.5,5,0.5), rep(c(0),each=length(seq(0.5,5,0.5))))
for (df in dfs) {
  cover <- poisson_expanding_window_sv(100,0.5,0.05,100,R=10,df=df)
  covers <- cbind(covers, (as.data.frame(cover[,-1])))
  print(covers)
}

# why is the cover getting worse as the window size increases?



