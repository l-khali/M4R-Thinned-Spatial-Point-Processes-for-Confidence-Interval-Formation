# REPEAT THIS WITH THE VARIANCE CALCULATED ANALYTICALLY!!!!


# # obtaining histogram of theta hat
W <- 1
A <- W^2
intensity <- 250
true_intensity <- intensity
intensity_hats <- c()
for (i in 1:10000) {
  p <- rpoispp(intensity,win=square(W))
  i <- npoints(p)/W^2
  intensity_hats <- c(intensity_hats, i)
}
hist(intensity_hats)
sigma2_spatial_intensity <- var((intensity_hats))
sigma2_spatial_intensity


# # theoretical var (values too large, returns infinity)
# A <- W^2
# n <- A*intensity
# sigma2_spatial_intensity <- (sqrt(n)/(thinning_param*A))^2 * (intensity*A)^n * exp(-intensity*A) * (thinning_param*(thinning_param+1)^(n-2)*(n*thinning_param+1))/(factorial(n-1)) - (intensity*A)^(2*n)*exp(-2*intensity*A)*(thinning_param^2*(thinning_param-1)^(2*(n-1)))/(factorial(n-1))^2


check_thinned_variance_spatial_intensity <- function(thinning_param, intensity) {
  thinned_intensity_hats <- c()
  normalised_thinned_intensity_hats <- c()
  for (i in 1:100) {
    p <- rpoispp(intensity,win=square(W))
    for (j in 1:500) {
      unif <- runif(npoints(p), 0, 1)
      thinned_p_df <- as.data.frame(p)[which(unif < thinning_param),]
      thinned_p <- as.ppp(thinned_p_df, W=square(W))
      i <- npoints(thinned_p)/(W^2*thinning_param)
      
      thinned_intensity_hats <- c(thinned_intensity_hats, i)
      normalised_thinned_intensity_hats <- c(normalised_thinned_intensity_hats, sqrt(npoints(p))*(i-true_intensity))
    }
  }
  sigma2thinned_intensity <- var(thinned_intensity_hats)
  return(sigma2thinned_intensity)
}

thinning_params <- seq(0.1,0.9,0.1)
thinned_sigmas_intensity <- c()
for (thinning_param in thinning_params) {
  print(thinning_param)
  thinned_sigmas_intensity <- c(thinned_sigmas_intensity, check_thinned_variance_spatial_intensity(thinning_param, intensity))
}

# par(mfrow = c(1, 2))
par(mfrow = c(1, 1))
# hopefully the lines are very close (not very good)
plot(thinning_params, thinned_sigmas_intensity,type="l",col=2,xlab=" ",ylab=" ",font.main=2,main=("Variance of Intensity Estimates"),lwd=1.5)
title(mgp=c(2.5,0,0),xlab=TeX("Thinning parameter, $p$"),ylab=TeX(r'(Sample Variance, $\hat{\sigma}$)'))
lines(thinning_params, sigma2_spatial_intensity/(thinning_params/(1-thinning_params)), col=3,lwd=1.5)
legend(0.7, 2450, c(TeX(r'($\hat{\sigma}/p$)'), TeX(r'($\hat{sigma}_s$)')), col = c(3, 2), lty = c(1, 1),
       merge = TRUE)
# lines(thinning_params, (true_k)*((1-thinning_params)/thinning_params))




