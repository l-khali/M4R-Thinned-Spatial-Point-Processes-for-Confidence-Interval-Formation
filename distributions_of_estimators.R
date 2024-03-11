# obtaining histogram of theta hat
n <-100000
theta_hats <- c()
for (i in 1:1000) {
  p <- rpois(n,10)
  theta_hats <- c(theta_hats, mean(p))
}
hist(theta_hats)
sigma2 <- var(sqrt(n)*(theta_hats-10))



thinned_theta_hats <- c()
normalised_thinned_theta_hats <- c()
thinning_param <- 0.5
p <- rpois(n,10)
for (i in 1:10000) {
  unif <- runif(length(p), 0, 1)
  thinned_p <- p[which(unif < thinning_param)]
  thinned_theta_hats <- c(thinned_theta_hats, mean(thinned_p))
  normalised_thinned_theta_hats <- c(normalised_thinned_theta_hats, sqrt(n)*(mean(thinned_p)-mean(rpois(n,10))))
}
hist(thinned_theta_hats)
# sigma2hat <- var(sqrt(n)*(thinned_theta_hats-mean(theta_hats)))
sigma2thinned <- var(normalised_thinned_theta_hats)
print(paste("sigma2hat",sigma2thinned))
print(paste("dsigma2/r",sigma2*((1-thinning_param)/thinning_param)))





# checking consistency
x <- 10
thinning_param <- 0.8
ns <- c(5,10,20,30,40)
thinned_hist <- c()
observed_hist <- c()
for (n in ns) {
  thinned <- c()
  observed <- c()
  for (sim in 1:1000) {
    for (i in 1:100) {
      p <- rpois(n, 10)
      observed <- c(observed, mean(p))
      unif <- runif(n, 0, 1)
      thinned_p <- p[which(unif < thinning_param)]
      thinned <- c(thinned, mean(thinned_p))
    }
    
    
    
    # thinned_prob <- length(thinned[which(thinned < x)])/length(thinned)
    # observed_prob <- length(observed[which(observed < x)])/length(observed)
    # thinned_hist <- c(thinned_hist, thinned_prob)
    # observed_hist <- c(observed_hist, observed_prob)
  }
  qqplot(thinned,observed)
}



plot((thinned_hist-observed_hist))



