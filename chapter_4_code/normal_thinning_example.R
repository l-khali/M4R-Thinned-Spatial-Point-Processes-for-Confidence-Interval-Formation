library(pals)
library(latex2exp)
library(spatstat)
palette(alphabet(26))

pois_thinning <- function(thinning_param = 0.5, N = 100, B = 1000, nsim=500) {
  
  bootstrapped_estimates <- c()
  thinned_estimates <- c()
  
  for (sim in 1:nsim) {
    population <- rpois(N,10)
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
  if (thinning_param < 0.5 ) {
    plot(thin_hist, col =rgb(116/255,10/255,255/255,alpha=0.7), freq=FALSE,xlab="",ylab="",main="",cex.lab=1.5,xlim=c(7,13),ylim=c(0,0.9))
    plot(boot_hist, col = rgb(43/255,206/255,72/255,alpha=0.7),freq=FALSE,add=TRUE)
    txt <- paste0("Poisson $lambda$ Estimates, $p = ",thinning_param,"$")
    title(TeX(txt), line=0.4, cex.main=1.5)
    # legend(-0.79, 3, legend=c("Bootstrapped Estimates", "Thinned Estimates"),
    #        col=c(rgb(43/255,206/255,72/255,alpha=0.7),rgb(116/255,10/255,255/255,alpha=0.7)), lty=1,lwd=3,cex=1.5)
  } else {
    plot(thin_hist, col = rgb(116/255,10/255,255/255,alpha=0.7), freq=FALSE,main="",xlab="",ylab="",cex.lab=1.5,xlim=c(8,12))
    plot(boot_hist, col = rgb(43/255,206/255,72/255,alpha=0.7), freq=FALSE, xlim=c(-0.7,0.7),main="", add=TRUE)
    txt <- paste0("Poisson $lambda$ Estimates, $p = ",thinning_param,"$")
    title(TeX(txt), line=0.4, cex.main=1.5)    }

  return(thinned_estimates)
}

# par(mfrow = c(1, 3))
# for (param in c(0.2,0.5,0.8)) {
#   normal_thinning(thinning_param = param, N=100)
# }

exp_thinning <- function(thinning_param = 0.5, N = 100, B = 1000, nsim=500) {
  
  bootstrapped_estimates <- c()
  thinned_estimates <- c()
  
  for (sim in 1:nsim) {
    population <- rexp(N, rate=10)
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
  if (thinning_param < 0.51 ) {
    plot(thin_hist, col =rgb(116/255,10/255,255/255,alpha=0.7), xlim=c(0.05,0.2), ylim=c(0,32),freq=FALSE,xlab="",ylab="",main="",cex.lab=1.5)
    plot(boot_hist, col = rgb(43/255,206/255,72/255,alpha=0.7),freq=FALSE,add=TRUE)
    txt <- paste0("Exponential $lambda$ Estimates, $p = ",thinning_param,"$")
    title(TeX(txt), line=0.4, cex.main=1.5)
    if (thinning_param < 0.5) {
      # legend(0.051, 32, legend=c("Bootstrapped Estimates", "Thinned Estimates"),
             # col=c(rgb(43/255,206/255,72/255,alpha=0.7),rgb(116/255,10/255,255/255,alpha=0.7)), lty=1,lwd=3,cex=1.3)
    }
  } else {
    plot(thin_hist, col = rgb(116/255,10/255,255/255,alpha=0.7),xlim=c(0.05,0.2), freq=FALSE,main="",cex.lab=1.5,xlab="",ylab="")
    plot(boot_hist, col = rgb(43/255,206/255,72/255,alpha=0.7), freq=FALSE,main="", add=TRUE)
    txt <- paste0("Exponential $lambda$ Estimates, $p = ",thinning_param,"$")
    title(TeX(txt), line=0.4, cex.main=1.5)    }
  
  return(thinned_estimates)
}
# 
# par(mfrow = c(1, 3))
# for (param in c(0.2,0.5,0.8)) {
#   exp_thinning(thinning_param = param, N=100)
# }



# combined plots
par(mfrow = c(2, 3))
op <- par(mfrow = c(2,3),
          oma = c(5,4,2,1) + 0.1,
          mar = c(2,2,2,2) + 0.1)

for (param in c(0.2,0.5,0.8)) {
  pois_thinning(thinning_param = param, N=100)
}

for (param in c(0.2,0.5,0.8)) {
  exp_thinning(thinning_param = param, N=100)
}

txt <- paste0("Histograms of Thinned and Bootstrapped Estimates")
title(TeX(txt),line=1,outer=TRUE, cex.main=1.6)
title(xlab="Estimates", ylab="Density", outer=TRUE, line=1, cex.lab=1.5)
legend(-0.2,90,legend = c("Thinned Estimates", "Bootstrapped"), col = c(rgb(116/255,10/255,255/255,alpha=0.7),rgb(43/255,206/255,72/255,alpha=0.7)), lwd = 5, horiz = TRUE, cex = 1.3, seg.len=1, bty = 'n',xpd=NA)
par(op)

