efron_2 <- function(thinning_param = 0.5, B = 50, n = 13) {
  population <- rnorm(n)
  theta_F <- median(population)
  sigma_F <- sqrt(var(population))
  # theta_F <- 0
  # sigma_F <- 1
  bootstrapped_Rs <- c()
  thinned_Rs <- c()
  
  for (b in 1:B) {
    bootstrap_sample <- sample(population, n, replace = TRUE)
    bootstrap_t <- median(bootstrap_sample)
    
    rand <- runif(n,0,1)
    thinned_sample <- population[which(rand < thinning_param)]
    thinned_t <- median(thinned_sample)
    
    # bootstrapped_Rs <- c(bootstrapped_Rs, abs(bootstrap_t - theta_F)/sigma_F)
    # thinned_Rs <- c(thinned_Rs, abs(thinned_t - theta_F)/sigma_F)
    bootstrapped_Rs <- c(bootstrapped_Rs, exp(bootstrap_t - theta_F))
    thinned_Rs <- c(thinned_Rs, exp(thinned_t - theta_F))
  }
  c1 <- rgb(255,0,0,max = 255, alpha = 90, names = "lt.blue")
  c2 <- rgb(0,255,0, max = 255, alpha = 90, names = "lt.pink")
  
  # boot_hist <- hist(bootstrapped_Rs, plot = FALSE)
  # thin_hist <- hist(thinned_Rs, plot = FALSE)
  # plot(boot_hist, col = c1, xlim=c(-0.5,0.5), ylim=c(0,15), freq=FALSE, main=paste("thinning parameter =", thinning_param))
  # plot(thin_hist, col = c2, add=TRUE, xlim=c(-0.5,0.5), ylim=c(0,20), freq=FALSE)
  qqplot(bootstrapped_Rs, thinned_Rs, main=paste("thinning parameter =", thinning_param))
  abline(0, 1, col = 'red')
}

par(mfrow = c(2, 5))
for (param in seq(0.1,1,0.1)) {
  efron_2(thinning_param = param, B=1000, n=50)
}