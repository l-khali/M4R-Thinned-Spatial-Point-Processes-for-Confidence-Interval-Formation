# obtaining histogram of theta hat
n <-100000
theta_hats <- c()
for (i in 1:1000) {
  p <- rpois(n,10)
  theta_hats <- c(theta_hats, mean(p))
}
hist(theta_hats)
sigma2 <- var(sqrt(n)*(theta_hats-10))


check_thinned_variance_1d <- function(thinning_param, intensity) {
  thinned_theta_hats <- c()
  normalised_thinned_theta_hats <- c()
  p <- rpois(n,intensity)
  for (i in 1:10000) {
    unif <- runif(length(p), 0, 1)
    thinned_p <- p[which(unif < thinning_param)]
    thinned_theta_hats <- c(thinned_theta_hats, mean(thinned_p))
    normalised_thinned_theta_hats <- c(normalised_thinned_theta_hats, sqrt(n)*(mean(thinned_p)-intensity))
  }
  sigma2thinned <- var(normalised_thinned_theta_hats)
  print(paste("sigma2hat",sigma2thinned))
  print(paste("dsigma2/r",sigma2*((1-thinning_param)/thinning_param)))
  return(sigma2thinned)
}

thinning_params <- seq(0.1,0.9,0.1)
thinned_sigmas_exp <- c()
for (thinning_param in thinning_params) {
  thinned_sigmas_exp <- c(thinned_sigmas_exp, check_thinned_variance_1d(thinning_param, 10))
}

# hopefully the lines are very close
plot(thinning_params, thinned_sigmas_exp,type="l")
lines(thinning_params, sigma2*((1-thinning_params)/thinning_params))


par(mfrow = c(1, 2))
plot(thinning_params, thinned_sigmas_exp,type="l",col=2,xlab="p",ylab=TeX("Thinned Variance, ${sigma}_s$"), main="Exponential Distribution")
lines(thinning_params, sigma2*((1-thinning_params)/thinning_params),col=3)
plot(thinning_params, thinned_sigmas,type="l",col=2,xlab="p",ylab=TeX("Thinned Variance, ${sigma}_s$"), main="Homogenous Poisson")
lines(thinning_params, sigma2_spatial*((1-thinning_params)/thinning_params),col=3)


















# # rubbish:
# 
# # checking consistency
# x <- 10
# thinning_param <- 0.8
# ns <- c(5,10,20,30,40)
# thinned_hist <- c()
# observed_hist <- c()
# for (n in ns) {
#   thinned <- c()
#   observed <- c()
#   for (sim in 1:1000) {
#     for (i in 1:100) {
#       p <- rpois(n, 10)
#       observed <- c(observed, mean(p))
#       unif <- runif(n, 0, 1)
#       thinned_p <- p[which(unif < thinning_param)]
#       thinned <- c(thinned, mean(thinned_p))
#     }
#     
#     
#     
#     # thinned_prob <- length(thinned[which(thinned < x)])/length(thinned)
#     # observed_prob <- length(observed[which(observed < x)])/length(observed)
#     # thinned_hist <- c(thinned_hist, thinned_prob)
#     # observed_hist <- c(observed_hist, observed_prob)
#   }
#   qqplot(thinned,observed)
# }
# 
# 
# 
# plot((thinned_hist-observed_hist))
# 
# 
# 
