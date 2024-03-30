library(spatstat)

# generate an inhomogenous point process with a simple intensity function
# intensity <- function(x,y) {
#   return(exp(1 + 3*x))
# }
# 
# p <- rpoispp(intensity, win=square(1))
# plot(p)

# attempting to estimate the coefficients of the intensity using the thinned sets
inhom_pois_parametric_intensity <- function(data, thinning_param, alpha, W=1, a, b, R=100) {
  
  # confidence intervals calculated using quantiles of samples
  processdf <- as.data.frame(data)
  # r <- seq(0.0, 0.14, 0.01)
  as <- c()
  bs <- c()
  for (i in 1:R) {
    unif <- runif(npoints(data), 0, 1)
    subprocessdf <- processdf[which(unif < thinning_param),]
    subp <- as.ppp(subprocessdf, W=square(W))
    fit <- ppm(subp ~ x)
    print(coef(fit))
    as <- c(as, coef(fit)[1]) # need to divide the intercept by thinning parameter?
    bs <- c(bs, coef(fit)[2])
  }
  
  fit <- ppm(data ~ x)
  ahat <- coef(fit)[1]
  bhat <- coef(fit)[2]
  scalar <- thinning_param
  t <- qt((1-alpha/2), df = R-1)
  
  a_lower_approx <- c()
  a_upper_approx <- c()
  
  b_lower_approx <- c()
  b_upper_approx <- c()
  
  a_lower_approx <- c(a_lower_approx, ahat - t*sqrt(var(as)*scalar))
  a_upper_approx <- c(a_upper_approx, ahat + t*sqrt(var(as)*scalar))
  
  b_lower_approx <- c(b_lower_approx, bhat - t*sqrt(var(bs)*scalar))
  b_upper_approx <- c(b_upper_approx, bhat + t*sqrt(var(bs)*scalar))
  
  lowers <- c(a_lower_approx, b_lower_approx)
  uppers <- c(a_upper_approx, b_upper_approx)
  
  return(cbind(lowers, uppers))
}

inhom_pois_parametric_thinning <- function(nsim, thinning_param, alpha, a, b, R=100) {
  
  intensity1 <- function(x,y) {
    return(exp(a+b*x))
  }
  
  Ws <- seq(2,5,0.5)
  cover <- rep(c(0),each=length(Ws))
  acoverage <- cbind(Ws, cover)
  bcoverage <- cbind(Ws, cover)
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity1, win=square(Ws[i]))
      confidences <- inhom_pois_parametric_intensity(p, thinning_param, alpha, W=Ws[i], a, b, R=R)
      # print(confidences)
      a_lower <- confidences[1,1]
      a_upper <- confidences[1,2]
      b_lower <- confidences[2,1]
      b_upper <- confidences[2,2]
      if (a_lower <= a & a <= a_upper) {
        acoverage[i,2] <- acoverage[i,2] + 1/nsim
      }
      if (b_lower <= b & b <= b_upper) {
        bcoverage[i,2] <- bcoverage[i,2] + 1/nsim
      }
    }
    print(acoverage[i,])
    print(bcoverage[i,])
  }
  return(cbind(acoverage,bcoverage[,2]))
  
}

parametric9 <- inhom_pois_parametric_thinning(10,0.9,0.05,1,3,R=50)

parametric5 <- inhom_pois_parametric_thinning(100,0.5,0.05,0.5,2,R=500)

parametric2 <- inhom_pois_parametric_thinning(10,0.2,0.05,1,3,R=50)
