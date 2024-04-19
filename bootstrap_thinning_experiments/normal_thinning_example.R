normal_thinning <- function(thinning_param = 0.5, N = 100, B = 1000, nsim=50) {
  
  bootstrapped_estimates <- c()
  thinned_estimates <- c()
  
  for (sim in 1:nsim) {
    population <- rnorm(N)
    for (b in 1:B) {
      s <- sample(population, N, replace = TRUE)
      bootstrapped_estimates <- c(bootstrapped_estimates, mean(s))
      
      rand <- runif(N,0,1)
      thinned <- population[which(rand < thinning_param)]
      thinned_estimates <- c(thinned_estimates, mean(thinned))
    }
  }
  # 
  c1 <- rgb(255,0,0,max = 255, alpha = 90, names = "lt.blue")
  c2 <- rgb(0,255,0, max = 255, alpha = 90, names = "lt.pink")

  # qqplot(bootstrapped_estimates,thinned_estimates,add.line=TRUE)
  # abline(lm(sort(thinned_estimates) ~ sort(bootstrapped_estimates)), col = "red", lwd = 2, lty = 2)
  boot_hist <- hist(bootstrapped_estimates, plot = FALSE)
  thin_hist <- hist(thinned_estimates, plot = FALSE)
  if (thinning_param < param ) {
    plot(boot_hist, col = 22, add=TRUE,freq=FALSE, xlim=c(-0.7,0.7),main=paste("N =", N, " nsim =", B, " thinning parameter =", thinning_param))
    plot(thin_hist, col = 7, xlim=c(-0.7,0.7), freq=FALSE)
  } else {
    plot(thin_hist, col = 7, add=TRUE,xlim=c(-0.7,0.7), freq=FALSE)
    plot(boot_hist, col = 22, freq=FALSE, xlim=c(-0.7,0.7),main=paste("N =", N, " nsim =", B, " thinning parameter =", thinning_param))
  }

  return(thinned_estimates)
}

palette(alphabet(26))

par(mfrow = c(1, 3))
for (param in c(0.2,0.5,0.8)) {
  normal_thinning(thinning_param = param, N=100)
}

legend(-1, 6, legend=c("bootstrapped", "thinned"),
       col=c(22,7), lty=1)

