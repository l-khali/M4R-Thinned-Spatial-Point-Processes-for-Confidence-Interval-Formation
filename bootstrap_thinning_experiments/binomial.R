binomial_experiment <- function(thinning_param = 0.5, N = 100, B = 10000, n=10) {
  population <- rbinom(N, n, 0.6)
  bootstrapped_estimates <- c()
  thinned_estimates <- c()
  for (b in 1:B) {
    s <- sample(population, N, replace = TRUE)
    bootstrapped_estimates <- c(bootstrapped_estimates, mean(s)/n)
    
    rand <- runif(N,0,1)
    thinned <- population[which(rand < thinning_param)]
    thinned_estimates <- c(thinned_estimates, mean(thinned)/n)
  }
  
  qqplot(bootstrapped_estimates, thinned_estimates, main=paste("Binomial, thinning parameter =", thinning_param))
  abline(0, 1, col = 'red')
}

par(mfrow = c(2, 5))
for (param in seq(0.1,1,0.1)) {
  binomial_experiment(thinning_param = param, N=100)
}
