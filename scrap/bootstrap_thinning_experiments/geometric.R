geometric_experiment <- function(thinning_param = 0.5, N = 100, B = 10000, prob=0.6) {
  population <- rgeom(N, prob=prob)
  bootstrapped_estimates <- c()
  thinned_estimates <- c()
  for (b in 1:B) {
    s <- sample(population, N, replace = TRUE)
    bootstrapped_estimates <- c(bootstrapped_estimates, 1/mean(s))
    
    rand <- runif(N,0,1)
    thinned <- population[which(rand < thinning_param)]
    thinned_estimates <- c(thinned_estimates, 1/mean(thinned))
  }
  
  qqplot(bootstrapped_estimates, thinned_estimates, main=paste("Geometric, thinning parameter =", thinning_param))
  abline(0, 1, col = 'red')
}

par(mfrow = c(2, 5))
for (param in seq(0.1,1,0.1)) {
  geometric_experiment(thinning_param = param, N=100)
}
