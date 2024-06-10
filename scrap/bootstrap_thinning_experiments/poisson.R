poisson_experiment <- function(thinning_param = 0.5, N = 100, B = 10000, lambda = 3) {
  population <- rpois(N, lambda=lambda)
  bootstrapped_estimates <- c()
  thinned_estimates <- c()
  for (b in 1:B) {
    s <- sample(population, N, replace = TRUE)
    bootstrapped_estimates <- c(bootstrapped_estimates, mean(s))
    
    rand <- runif(N,0,1)
    thinned <- population[which(rand < thinning_param)]
    thinned_estimates <- c(thinned_estimates, mean(thinned))
  }
  return(thinned_estimates)
  # qqplot(bootstrapped_estimates, thinned_estimates, main=paste("Poisson, thinning parameter =", thinning_param))
  # abline(0, 1, col = 'red')
  
}

par(mfrow = c(2, 5))
for (param in seq(0.1,1,0.1)) {
  poisson_experiment(thinning_param = param, N=100)
}
