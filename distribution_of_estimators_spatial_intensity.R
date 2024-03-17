# REPEAT THIS WITH THE VARIANCE CALCULATED ANALYTICALLY!!!!


# # obtaining histogram of theta hat
W <- 2
A <- W^2
intensity <- 50
true_intensity <- intensity
# intensity_hats <- c()
# for (i in 1:10000) {
#   p <- rpoispp(intensity,win=square(W))
#   i <- npoints(p)/W^2
#   intensity_hats <- c(intensity_hats, i)
# }
# hist(intensity_hats)
# sigma2_spatial_intensity <- var(sqrt(npoints(p))*(intensity_hats-true_intensity))
# sigma2_spatial_intensity


# theoretical var (values too large, returns infinity)
A <- W^2
n <- A*intensity
sigma2_spatial_intensity <- (sqrt(n)/(thinning_param*A))^2 * (intensity*A)^n * exp(-intensity*A) * (thinning_param*(thinning_param+1)^(n-2)*(n*thinning_param+1))/(factorial(n-1)) - (intensity*A)^(2*n)*exp(-2*intensity*A)*(thinning_param^2*(thinning_param-1)^(2*(n-1)))/(factorial(n-1))^2


check_thinned_variance_spatial_intensity <- function(thinning_param, intensity) {
  thinned_intensity_hats <- c()
  normalised_thinned_intensity_hats <- c()
  for (i in 1:200) {
    p <- rpoispp(intensity,win=square(W))
    for (j in 1:50) {
      unif <- runif(npoints(p), 0, 1)
      thinned_p_df <- as.data.frame(p)[which(unif < thinning_param),]
      thinned_p <- as.ppp(thinned_p_df, W=square(W))
      i <- npoints(thinned_p)/(W^2*thinning_param)
      
      thinned_intensity_hats <- c(thinned_intensity_hats, i)
      normalised_thinned_intensity_hats <- c(normalised_thinned_intensity_hats, sqrt(npoints(p))*(i-true_intensity))
    }
  }
  hist(thinned_intensity_hats)
  sigma2thinned_intensity <- var(normalised_thinned_intensity_hats)
  print(paste("sigma2hat",sigma2thinned_intensity))
  print(paste("dsigma2/r",sigma2_spatial_intensity*((1-thinning_param)/thinning_param)))
  return(sigma2thinned_intensity)
}

thinning_params <- seq(0.1,0.9,0.1)
thinned_sigmas_intensity <- c()
for (thinning_param in thinning_params) {
  thinned_sigmas_intensity <- c(thinned_sigmas_intensity, check_thinned_variance_spatial_intensity(thinning_param, intensity))
}

par(mfrow = c(1, 2))
par(mfrow = c(1, 1))
# hopefully the lines are very close (not very good)
plot(thinning_params, thinned_sigmas_intensity,type="l")
lines(thinning_params, sigma2_spatial_intensity*((1-thinning_params)/thinning_params))
# lines(thinning_params, (true_k)*((1-thinning_params)/thinning_params))

