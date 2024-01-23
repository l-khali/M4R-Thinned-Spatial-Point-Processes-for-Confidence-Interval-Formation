normal_thinning <- function(thinning_param = 0.5, N = 100, B = 1000) {
  population <- rnorm(N)
  bootstrapped_estimates <- c()
  thinned_estimates <- c()
  for (b in 1:B) {
    s <- sample(population, N, replace = TRUE)
    bootstrapped_estimates <- c(bootstrapped_estimates, mean(s))
    
    rand <- runif(N,0,1)
    thinned <- population[which(rand < thinning_param)]
    thinned_estimates <- c(thinned_estimates, mean(thinned))
  }
  # 
  # c1 <- rgb(255,0,0,max = 255, alpha = 90, names = "lt.blue")
  # c2 <- rgb(0,255,0, max = 255, alpha = 90, names = "lt.pink")
  # 
  # boot_hist <- hist(bootstrapped_estimates, plot = FALSE)
  # thin_hist <- hist(thinned_estimates, plot = FALSE)
  # plot(boot_hist, col = c1, xlim=c(-0.5,0.5), ylim=c(0,15), freq=FALSE, main=paste("N =", N, " nsim =", B, " thinning parameter =", thinning_param))
  # plot(thin_hist, col = c2, add=TRUE, xlim=c(-0.5,0.5), ylim=c(0,20), freq=FALSE)
  # legend(-1, 6, legend=c("bootstrapped", "thinned"),
         # col=c(c1,c2), lty=1)
  return(thinned_estimates)
}

par(mfrow = c(2, 5))
for (param in seq(0.1,1,0.1)) {
  normal_thinning(thinning_param = param, N=1000)
}
