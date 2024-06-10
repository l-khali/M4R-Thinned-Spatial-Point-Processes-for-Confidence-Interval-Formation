# obtaining histogram of theta hat
r <- c(0.0,0.05)
W <- 2
parent_intensity <- 25
matern_radius <- 0.1
daughter_intensity <- 10
# true_k <- pi*r[2]*r[2]
k_hats <- c()
for (i in 1:10000) {
  p <- rMatClust(parent_intensity,matern_radius,daughter_intensity,win=square(W))
  k <- Kest(p, r=r, correction=c("isotropic"))
  k_hats <- c(k_hats, k$iso[2])
}
hist(k_hats)
sigma2_spatial <- var((k_hats))
sigma2_spatial


check_thinned_variance_spatial <- function(thinning_param, intensity) {
  thinned_k_hats <- c()
  normalised_thinned_k_hats <- c()
  for (i in 1:200) {
    p <- rMatClust(parent_intensity,matern_radius,daughter_intensity,win=square(W))
    for (j in 1:50) {
      unif <- runif(npoints(p), 0, 1)
      thinned_p_df <- as.data.frame(p)[which(unif < thinning_param),]
      thinned_p <- as.ppp(thinned_p_df, W=square(W))
      k <- Kest(thinned_p, r=r, correction=c("isotropic"))
      thinned_k_hats <- c(thinned_k_hats, k$iso[2])
      normalised_thinned_k_hats <- c(normalised_thinned_k_hats, sqrt(npoints(p))*(k$iso[2]-true_k))
    }
  }
  sigma2thinned <- var(thinned_k_hats)
  print(paste("sigma2hat",sigma2thinned))
  print(paste("dsigma2/r",sigma2_spatial*((1-thinning_param)/thinning_param)))
  return(sigma2thinned)
}

thinning_params <- seq(0.1,0.9,0.1)
thinned_sigmas <- c()
for (thinning_param in thinning_params) {
  thinned_sigmas <- c(thinned_sigmas, check_thinned_variance_spatial(thinning_param, intensity))
}

par(mfrow = c(1, 1))
# hopefully the lines are very close (not very good)
plot(thinning_params, 250*thinned_sigmas,type="l")
lines(thinning_params, 250*sigma2_spatial/thinning_params)

scalars <- sigma2_spatial/thinned_sigmas
