library(spatstat)

# attempting to estimate the coefficients of the intensity using the thinned sets
inhom_pois_parametric_intensity <- function(data, thinning_param, alpha, W=1, R=100) {
  
  # confidence intervals calculated using quantiles of samples
  # processdf <- as.data.frame(data)
  as <- c()
  bs <- c()
  for (j in 1:R) {
    # unif <- runif(npoints(data), 0, 1)
    # subprocessdf <- processdf[which(unif < thinning_param),]
    # subp <- as.ppp(subprocessdf, W=square(W))
    subp <- rthin(data, thinning_param)
    fit <- ppm(subp ~ x)
    as <- c(as, coef(fit)[1]) # need to divide the intercept by thinning parameter?
    bs <- c(bs, coef(fit)[2])
  }
  
  fit <- ppm(data ~ x)
  ahat <- coef(fit)[1]
  bhat <- coef(fit)[2]
  scalar <- thinning_param/(1-thinning_param)
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
  
  Ws <- seq(0.5,2,0.25)
  cover <- rep(c(0),each=length(Ws))
  acoverage <- cbind(Ws, cover)
  bcoverage <- cbind(Ws, cover)
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity1, win=square(Ws[i]))
      confidences <- inhom_pois_parametric_intensity(p, thinning_param, alpha, W=Ws[i], R=R)
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

parametric8 <- inhom_pois_parametric_thinning(100,0.8,0.05,3,3,R=500)

parametric5 <- inhom_pois_parametric_thinning(100,0.5,0.05,3,3,R=500)

parametric2 <- inhom_pois_parametric_thinning(100,0.2,0.05,3,3,R=500)

save(parametric8,file="covers_ppm_intensity_8.RData")
save(parametric5,file="covers_ppm_intensity_5.RData")
save(parametric2,file="covers_ppm_intensity_2.RData")

par(mfrow=c(1,2))

plot(seq(0.5,2,0.25),parametric8[,2],type="l",col=7,lwd=1.5,xlab=TeX('Window Side Length, $W$'),ylab="Confidence Interval Cover",ylim=c(0.8,1))
lines(seq(0.5,2,0.25),parametric5[,2],ylim=c(0.5,1),col=15,lwd=1.5)
lines(seq(0.5,2,0.25),parametric2[,2],ylim=c(0.5,1),col=22,lwd=1.5)
abline(h=0.95, col=18,lwd=1.5,lty=2)
title(TeX("Parametric Est. of $a$"),line=0.8)

plot(seq(0.5,2,0.25),parametric8[,3],type="l",col=7,lwd=1.5,xlab=TeX('Window Side Length, $W$'),ylab="Confidence Interval Cover",ylim=c(0.8,1))
lines(seq(0.5,2,0.25),parametric5[,3],ylim=c(0.5,1),col=15,lwd=1.5)
lines(seq(0.5,2,0.25),parametric2[,3],ylim=c(0.5,1),col=22,lwd=1.5)
abline(h=0.95, col=18,lwd=1.5,lty=2)
legend(1.4,0.85,c("p=0.8","p=0.5","p=0.2"),col=c(7,15,22),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),cex=0.9)
title(TeX("Parametric Est. of $b$"),line=0.8)

title("Studentized CI Cover: Parametric Intensity of Spatial Point Processes", line = -2, outer = TRUE,cex.main=1.2)



