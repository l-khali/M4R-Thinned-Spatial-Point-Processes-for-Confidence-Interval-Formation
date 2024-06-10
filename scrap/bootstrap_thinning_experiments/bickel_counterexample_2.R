library(spatstat)

counterexample2 <- function(thinning_param = 0.5, N = 1000, B = 1000) {
  ns <- seq(N/10,N,N/10)
  bootstrapped_p <- c()
  thinned_p <- c()
  for (n in ns) {
    population <- runif(n)
    Xn <- max(population)
    bootstrapped_estimates <- c()
    thinned_estimates <- c()
    for (b in 1:B) {
      s <- sample(population, n, replace = TRUE)
      bootstrapped_estimates <- c(bootstrapped_estimates, n * (Xn - max(s)))
      
      rand <- runif(n,0,1)
      thinned <- population[which(rand < thinning_param)]
      thinned_estimates <- c(thinned_estimates, n * (Xn - max(thinned)))
    }
    bootstrapped_p <- c(bootstrapped_p,length(which(bootstrapped_estimates == 0)) / B)
    thinned_p <- c(thinned_p,length(which(thinned_estimates == 0)) / B)
  }
  c1 <- rgb(255,0,0,max = 255, alpha = 90, names = "lt.blue")
  c2 <- rgb(0,255,0, max = 255, alpha = 90, names = "lt.pink")
  
  plot(ns, bootstrapped_p, ylim=c(0,1),type="l", col=2, main=paste("Estimating maximum of uniform RV, thinning_param=", thinning_param))
  lines(ns, thinned_p, col=3)
  legend(2000, 1, col=c(2,3), legend=c("bootstrap","thinned"), lty=1)
  
  print(bootstrapped_p)
  print(thinned_p)
}

par(mfrow = c(2, 5))
for (param in seq(0.1,1,0.1)) {
  counterexample2(thinning_param = param, N=10000)
}
