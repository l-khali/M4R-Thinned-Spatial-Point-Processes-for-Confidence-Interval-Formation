library(spatstat)

# attempting to estimate the coefficients of the intensity using the thinned sets
inhom_pois_parametric_intensity_clustered <- function(data, thinning_param, alpha, W=1, R=100) {
  
  # confidence intervals calculated using quantiles of samples
  # processdf <- as.data.frame(data)
  kappas <- c()
  scales <- c()
  for (j in 1:R) {
    subp <- rthin(data, thinning_param)
    fit <- kppm(subp ~ x, "Thomas", method="clik2")
    kappas <- c(kappas, fit$par[1]) # need to divide the intercept by thinning parameter?
    scales <- c(scales, fit$par[2])
  }
  
  fit <- kppm(data ~ x, "Thomas", method="clik2")
  kappahat <- fit$par[1]
  scalehat <- fit$par[2]
  scalar <- thinning_param/(1-thinning_param)
  t <- qt((1-alpha/2), df = R-1)
  
  kappa_lower_approx <- c()
  kappa_upper_approx <- c()
  
  scale_lower_approx <- c()
  scale_upper_approx <- c()
  
  kappa_lower_approx <- c(kappa_lower_approx, kappahat - t*sqrt(var(kappas)*scalar))
  kappa_upper_approx <- c(kappa_upper_approx, kappahat + t*sqrt(var(kappas)*scalar))
  
  scale_lower_approx <- c(scale_lower_approx, scalehat - t*sqrt(var(scales)*scalar))
  scale_upper_approx <- c(scale_upper_approx, scalehat + t*sqrt(var(scales)*scalar))
  
  lowers <- c(kappa_lower_approx, scale_lower_approx)
  uppers <- c(kappa_upper_approx, scale_upper_approx)
  
  return(cbind(lowers, uppers))
}

inhom_pois_parametric_thinning_clustered <- function(nsim, thinning_param, alpha, kappa, scale, mu, R=100) {
  
  Ws <- seq(3,5.5,1)
  cover <- rep(c(0),each=length(Ws))
  kappacoverage <- cbind(Ws, cover)
  scalecoverage <- cbind(Ws, cover)
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rThomas(kappa, scale, mu, win=square(Ws[i]))
      confidences <- inhom_pois_parametric_intensity_clustered(p, thinning_param, alpha, W=Ws[i], R=R)
      kappa_lower <- confidences[1,1]
      kappa_upper <- confidences[1,2]
      scale_lower <- confidences[2,1]
      scale_upper <- confidences[2,2]
      if (kappa_lower <= kappa & kappa <= kappa_upper) {
        kappacoverage[i,2] <- kappacoverage[i,2] + 1/nsim
      }
      if (scale_lower <= scale^2 & scale^2 <= scale_upper) {
        scalecoverage[i,2] <- scalecoverage[i,2] + 1/nsim
      }
    }
    print(kappacoverage[i,])
    print(scalecoverage[i,])
  }
  return(cbind(kappacoverage,scalecoverage[,2]))
  
}

parametric8_thomas <- inhom_pois_parametric_thinning_clustered(100,0.8,0.05,25,0.1,10, R=500)
save(parametric8_thomas,file="covers_ppm_intensity_8_thomas")
parametric5_thomas <- inhom_pois_parametric_thinning_clustered(100,0.5,0.05,25,0.1,10, R=500)
save(parametric5_thomas,file="covers_ppm_intensity_5_thomas")
parametric2_thomas <- inhom_pois_parametric_thinning_clustered(100,0.2,0.05,25,0.1,10, R=500)
save(parametric2_thomas,file="covers_ppm_intensity_2_thomas")

parametric8_thomas_cont <- inhom_pois_parametric_thinning_clustered(100,0.8,0.05,25,0.1,10, R=500)
save(parametric8_thomas_cont,file="covers_ppm_intensity_8_thomas_cont")
