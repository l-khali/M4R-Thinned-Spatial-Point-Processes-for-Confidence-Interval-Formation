library(spatstat)

# attempting to estimate the coefficients of the intensity using the thinned sets
inhom_pois_parametric_intensity_poly <- function(data, thinning_param, alpha, R=100) {
  
  # confidence intervals calculated using quantiles of samples
  # processdf <- as.data.frame(data)
  as <- c()
  bs <- c()
  cs <- c()
  # ds <- c()
  for (j in 1:R) {
    # unif <- runif(npoints(data), 0, 1)
    # subprocessdf <- processdf[which(unif < thinning_param),]
    # subp <- as.ppp(subprocessdf, W=square(W))
    subp <- rthin(data, thinning_param)
    fit <- ppm(subp ~ polynom(x,2))
    as <- c(as, coef(fit)[1]) # need to divide the intercept by thinning parameter?
    bs <- c(bs, coef(fit)[2])
    cs <- c(cs, coef(fit)[3])
    # ds <- c(ds, coef(fit)[4])
  }
  
  fit <- ppm(data ~ polynom(x,3))
  ahat <- coef(fit)[1]
  bhat <- coef(fit)[2]
  chat <- coef(fit)[3]
  # dhat <- coef(fit)[4]
  scalar <- thinning_param/(1-thinning_param)
  t <- qt((1-alpha/2), df = R-1)
  
  a_lower_approx <- c()
  a_upper_approx <- c()
  
  b_lower_approx <- c()
  b_upper_approx <- c()
  
  c_lower_approx <- c()
  c_upper_approx <- c()
  
  # d_lower_approx <- c()
  # d_upper_approx <- c()
  
  a_lower_approx <- c(a_lower_approx, ahat - t*sqrt(var(as)*scalar))
  a_upper_approx <- c(a_upper_approx, ahat + t*sqrt(var(as)*scalar))
  
  b_lower_approx <- c(b_lower_approx, bhat - t*sqrt(var(bs)*scalar))
  b_upper_approx <- c(b_upper_approx, bhat + t*sqrt(var(bs)*scalar))
  
  c_lower_approx <- c(c_lower_approx, ahat - t*sqrt(var(cs)*scalar))
  c_upper_approx <- c(c_upper_approx, ahat + t*sqrt(var(cs)*scalar))
  
  # d_lower_approx <- c(b_lower_approx, bhat - t*sqrt(var(ds)*scalar))
  # d_upper_approx <- c(b_upper_approx, bhat + t*sqrt(var(ds)*scalar))
  
  lowers <- c(a_lower_approx, b_lower_approx, c_lower_approx)
  uppers <- c(a_upper_approx, b_upper_approx, c_upper_approx)
  
  return(cbind(lowers, uppers))
}

inhom_pois_parametric_thinning_poly <- function(nsim, thinning_param, alpha, a, b, c, R=100) {
  
  intensity1 <- function(x,y) {
    return(exp(a+b*x+c*x^2))
  }
  
  Ws <- seq(1,4,0.5)
  cover <- rep(c(0),each=length(Ws))
  acoverage <- cbind(Ws, cover)
  bcoverage <- cbind(Ws, cover)
  ccoverage <- cbind(Ws, cover)
  # dcoverage <- cbind(Ws, cover)
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity1, win=square(Ws[i]))
      confidences <- inhom_pois_parametric_intensity_poly(p, thinning_param, alpha, R=R)
      print(confidences)
      a_lower <- confidences[1,1]
      a_upper <- confidences[1,2]
      b_lower <- confidences[2,1]
      b_upper <- confidences[2,2]
      c_lower <- confidences[3,1]
      c_upper <- confidences[3,2]
      # d_lower <- confidences[4,1]
      # d_upper <- confidences[4,2]
      if (a_lower <= a & a <= a_upper) {
        acoverage[i,2] <- acoverage[i,2] + 1/nsim
      }
      if (b_lower <= b & b <= b_upper) {
        bcoverage[i,2] <- bcoverage[i,2] + 1/nsim
      }
      if (c_lower <= c & c <= c_upper) {
        ccoverage[i,2] <- ccoverage[i,2] + 1/nsim
      }
      # if (d_lower <= d & d <= d_upper) {
      #   dcoverage[i,2] <- dcoverage[i,2] + 1/nsim
      # }
    }
    print(acoverage[i,])
    print(bcoverage[i,])
    print(ccoverage[i,])
    # print(dcoverage[i,])
  }
  return(cbind(acoverage,bcoverage[,2],ccoverage[,2]))
  
}

parametric8_poly <- inhom_pois_parametric_thinning_poly(300,0.8,0.05,1,0.5,0.5,0.5,R=500)

parametric5_poly <- inhom_pois_parametric_thinning_poly(100,0.5,0.05,3,0.5,0.5,R=100)

parametric2_poly <- inhom_pois_parametric_thinning_poly(300,0.2,0.05,3,3,1,1,R=500)

save(parametric8,file="covers_ppm_intensity_8")
save(parametric5,file="covers_ppm_intensity_5")
save(parametric2,file="covers_ppm_intensity_2")
