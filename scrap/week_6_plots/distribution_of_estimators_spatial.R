purple <- rgb(116/255,10/255,255/255)
green <- rgb(43/255,206/255,72/255)

# obtaining histogram of theta hat
r <- c(0.0,0.1)
W <- 2
intensity <- 2500
true_k <- pi*r[2]*r[2]
k_hats <- c()
for (i in 1:10000) {
  p <- rpoispp(intensity,win=square(W))
  k <- Kest(p, r=r, correction=c("isotropic"))
  k_hats <- c(k_hats, k$iso[2])
}
hist(k_hats)
sigma2_spatial <- var(sqrt(intensity*W^2)*(k_hats-true_k))
sigma2_spatial


check_thinned_variance_spatial <- function(thinning_param, intensity) {
  thinned_k_hats <- c()
  normalised_thinned_k_hats <- c()
  for (i in 1:1) {
    p <- rpoispp(intensity,win=square(W))
    for (j in 1:1000) {
      k_hat <- Kest(p, r=r, correction=c("isotropic"))$iso[2]
      unif <- runif(npoints(p), 0, 1)
      thinned_p_df <- as.data.frame(p)[which(unif < thinning_param),]
      thinned_p <- as.ppp(thinned_p_df, W=square(W))
      k <- Kest(thinned_p, r=r, correction=c("isotropic"))
      thinned_k_hats <- c(thinned_k_hats, k$iso[2])
      normalised_thinned_k_hats <- c(normalised_thinned_k_hats, sqrt(npoints(p))*(k$iso[2]-k_hat))
    }
  }
  sigma2thinned <- var(normalised_thinned_k_hats)
  print(paste("sigma2hat",sigma2thinned))
  print(paste("dsigma2/r",sigma2_spatial*((1-thinning_param)/thinning_param)))
  return(sigma2thinned)
}

thinning_params <- seq(0.1,0.9,0.1)
thinned_sigmas <- c()
for (thinning_param in thinning_params) {
  thinned_sigmas <- c(thinned_sigmas, check_thinned_variance_spatial(thinning_param, intensity))
}

par(mfrow = c(1, 2))
# hopefully the lines are very close
plot(thinning_params, thinned_sigmas,type="l",xlab=" ",ylab=" ",font.main=2,main=("Variance of Ripley's K Estimates on \n Poisson Point Process"),lwd=2,col=purple)
title(mgp=c(2.5,0,0),xlab=TeX("Thinning parameter, $p$"),ylab=TeX(r'(Sample Variance, $\hat{\sigma}$)'))
lines(thinning_params, sigma2_spatial*(1-thinning_params)/thinning_params, lwd=2,col=green)
legend(0.55, 0.0025, c(TeX(r'($\hat{\sigma}_s$)'), TeX(r'($\hat{sigma}(1-p)/p$)')), col = c(purple, green), lwd = c(2, 2),
       merge = TRUE)
# abline(h=sigma2_spatial)
# lines(thinning_params, (true_k)*((1-thinning_params)/thinning_params))

# lagache scaling
plot(thinning_params, thinned_sigmas,type="l",xlab=" ",ylab=" ",font.main=2,main=("Variance of Ripley's K Estimates on \n Poisson Point Process"),lwd=2,col=purple)
title(mgp=c(2.5,0,0),xlab=TeX("Thinning parameter, $p$"),ylab=TeX(r'(Sample Variance, $\hat{\sigma}$)'))
lines(thinning_params, sigma2_spatial/thinning_params^2, lwd=2,col=green)
legend(0.55, 0.0025, c(TeX(r'($\hat{\sigma}_s$)'), TeX(r'($\hat{sigma}/p^2$)')), col = c(purple, green), lwd = c(2, 2),
       merge = TRUE)


# theoretic variance from lagache + lang paper (works well!)
area <- W^2
boundary <- 4*W
beta <- pi*r[2]^2/area
gamma <- boundary*r[2]/area
sigma2_theos <- c()
sigma2_theos_n <- c()

for (thinning_param in thinning_params) {
  n <- intensity * area
  np <- intensity*area*thinning_param
  sigma2_theo <- n*((2*area^2*beta)/(np^2))*(1+0.305*gamma+beta*(-1+0.0132*np*gamma))
  sigma2_theos <- c(sigma2_theos, sigma2_theo)
  
  sigma2_theo_n <- n*((2*area^2*beta)/(n^2))*(1+0.305*gamma+beta*(-1+0.0132*n*gamma))
  sigma2_theos_n <- c(sigma2_theos_n, sigma2_theo_n/(thinning_param^2))
}

# plot(thinning_params, thinned_sigmas,type="l", col=2, xlab=TeX("$p$"))
# lines(thinning_params, sigma2_theos,col=3)
# 
# 
# plot(thinning_params, thinned_sigmas,type="l", col=2, xlab=TeX("Thinning parameter, $p$"),ylab=TeX("Thinned Variance, ${sigma}_s$"),main=TeX("$Var(K_s)$ on Homogenous Poisson"))
# lines(thinning_params, sigma2_theos_n,col=3)
# 
# plot(thinning_params, sigma2_theos,type="l", col=4, xlab=TeX("$p$"))
# lines(thinning_params, sigma2_theos_n,col=3)


plot(thinning_params, thinned_sigmas,type="l", col=2, xlab=TeX("Thinning parameter, $p$"),ylab=TeX("Thinned Variance, ${sigma}_s$"),main=TeX("$Var(K_s)$ on Homogenous Poisson"))
lines(thinning_params, sigma2_spatial*(1-thinning_params)/(thinning_params),col=3)

