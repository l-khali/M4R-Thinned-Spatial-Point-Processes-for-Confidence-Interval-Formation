# use histograms to show the distribution of normalised thinned estimates 
# is standard normal

# Thinning function for one dimensional process which statistic being intensity
thinning_thm3 <- function(data, thinning_param, alpha, R=99) {
  lambdas <- c()
  means <- c()
  for (i in 1:R) {
    # thinning the given process
    unif <- runif(length(data), 0, 1)
    subprocess_df <- data[which(unif < thinning_param)]
    r <- length(subprocess_df)
    d <- length(data) - r
    # calculating jackknife histogram
    j <- (mean(subprocess_df) - mean(data))*sqrt(length(data) * r/d)
    means <- c(means,mean(subprocess_df))
    lambdas <- c(lambdas, j)
  }
  lambdas[is.nan(lambdas)] <- 0 # subprocesses can have 0 points with low thinning param
  means[is.nan(means)] <- 0
  print(means)
  lambdas <- lambdas / sqrt(50)
  return(lambdas)
}

# grid of plots
par(mfrow = c(3, 3))
op <- par(mfrow = c(3,3),
          oma = c(5,4,2,1) + 0.1,
          mar = c(2,2,2,2) + 0.1)

# Homogenous Poisson Process (independent)
for (thinning_param in c(0.2,0.5,0.8)) {
  R <- 2000
  ns <- c(10,20,500)
  for (n in ns) {
    # homogenous poisson with mean 50
    p <- rpois(n,50)
    npoint <- nrow(as.data.frame(p))
    # obtaining jackknife histogram
    thinned_estimates <- thinning_thm3(p, thinning_param, 0.05, R=R)
    normal_dist <- rnorm(R, 0, 1)
    # histograms to observe convergence in distribution
    standard <- hist(normal_dist, plot = FALSE)
    thin_hist <- hist(thinned_estimates, plot = FALSE)
    if (n<1000) {
      plot(thin_hist, col =rgb(116/255,10/255,255/255,alpha=0.7), freq=FALSE,main="",xlab="",ylab="",cex.lab=1.5,ylim=c(0,0.55),xlim=c(-3,3))
      plot(standard, col = rgb(43/255,206/255,72/255,alpha=0.7),freq=FALSE,main="", ylab = "", xlab="",add=TRUE)
      txt <- paste0("$n = ",n,", p = ",thinning_param,"$")
      title(TeX(txt), line=0.5, cex.main=1.5)
    } else{
      plot(thin_hist, col =rgb(116/255,10/255,255/255,alpha=0.7), freq=FALSE,main="",xlab="",ylab="",cex.lab=1.5)
      plot(standard, col = rgb(43/255,206/255,72/255,alpha=0.7),freq=FALSE,main="", ylab = "", xlab="",add=TRUE)
      txt <- paste0("$n = ",n,", p = ",thinning_param,"$")
      title(TeX(txt), line=0.5, cex.main=1.5)
    }

  }
}
txt <- paste0("Estimating $lambda$ of Poisson Distribution")
title(TeX(txt),line=1,outer=TRUE, cex.main=1.6)
title(xlab="Estimates", ylab="Density", outer=TRUE, line=1, cex.lab=1.5)
legend(-15.6,2.38,legend = c("Thinned Estimates", "Standard Normal"), col = c(rgb(116/255,10/255,255/255,alpha=0.7),rgb(43/255,206/255,72/255,alpha=0.7)), lwd = 5, horiz = TRUE, cex = 1.3, seg.len=1, bty = 'n',xpd=NA)
par(op)

