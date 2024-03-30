# obtaining histogram of theta hat
r <- c(0.0,0.05)
W <- 2
parent_intensity <- 25
thomas_radius <- 0.1
daughter_intensity <- 10
true_k <- pi*r[2]*r[2] + (1/parent_intensity)*(1-exp((-r[2]^2)/(4*thomas_radius^2)))
k_hats <- c()
for (i in 1:10000) {
  p <- rThomas(parent_intensity,thomas_radius,daughter_intensity,win=square(W))
  k <- Kest(p, r=r, correction=c("isotropic"))
  k_hats <- c(k_hats, k$iso[2])
}
hist(k_hats)
sigma2_spatial_thomas <- var(sqrt(250)*(k_hats-true_k))
sigma2_spatial_thomas


check_thinned_variance_spatial_thomas <- function(thinning_param, intensity) {
  thinned_k_hats <- c()
  normalised_thinned_k_hats <- c()
  for (i in 1:200) {
    p <- rThomas(parent_intensity,thomas_radius,daughter_intensity,win=square(W))
    for (j in 1:50) {
      unif <- runif(npoints(p), 0, 1)
      thinned_p_df <- as.data.frame(p)[which(unif < thinning_param),]
      thinned_p <- as.ppp(thinned_p_df, W=square(W))
      k <- Kest(thinned_p, r=r, correction=c("isotropic"))
      thinned_k_hats <- c(thinned_k_hats, k$iso[2])
      normalised_thinned_k_hats <- c(normalised_thinned_k_hats, sqrt(npoints(p))*(k$iso[2]-true_k))
    }
  }
  sigma2thinned <- var(normalised_thinned_k_hats)
  print(paste("sigma2hat",sigma2thinned))
  print(paste("dsigma2/r",sigma2_spatial*((1-thinning_param)/thinning_param)))
  return(sigma2thinned)
}

thinning_params <- seq(0.1,0.9,0.1)
thinned_sigmas_thomas <- c()
for (thinning_param in thinning_params) {
  thinned_sigmas_thomas <- c(thinned_sigmas_thomas, check_thinned_variance_spatial_thomas(thinning_param, intensity))
}

plot(thinning_params,thinned_sigmas_thomas,type="l",col=2,lwd=1.5,main="Variance of Ripleys K on Thomas Process")
lines(thinning_params,sigma2_spatial_thomas/(thinning_params*(1-thinning_params)),lwd=1.5)

scalers <- sigma2_spatial_thomas/thinned_sigmas_thomas

