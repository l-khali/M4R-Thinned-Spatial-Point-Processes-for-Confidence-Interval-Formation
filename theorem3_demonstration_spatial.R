# Thinning function for homogenous poisson pp which statistic being the intensity
thinning_thm3_spatial <- function(data, thinning_param, alpha, R=99, inten, w) {
  Ks <- c()
  for (i in 1:R) {
    # thinning the given process
    subprocess_df <- rthin(data, thinning_param)
    r <- npoints(subprocess_df)
    d <- npoints(data) - r
    # calculating histogram
    j <- ((r/(thinning_param*w^2)) - (npoints(data)/(w^2)))*sqrt(npoints(data) * r/d)
    Ks <- c(Ks, j)
  }
  Ks[is.nan(Ks)] <- 0 # subprocesses can have 0 points with low thinning param
  Ks <- Ks / sqrt(inten/(w^2))
  return(Ks)
}

# grid of plots
par(mfrow = c(3, 3))
op <- par(mfrow = c(3,3),
          oma = c(5,4,2,1) + 0.1,
          mar = c(2,2,2,2) + 0.1)

# Homogenous Poisson Process (independent)
for (thinning_param in c(0.2,0.5,0.8)) {
  R <- 2000
  ws <- c(0.05,2,5)
  for (w in ws) {
    # homogenous poisson with mean 50
    p <- rpoispp(250,win=square(w))
    npoint <- nrow(as.data.frame(p))
    # obtaining jackknife histogram
    thinned_estimates <- thinning_thm3_spatial(p, thinning_param, 0.05, R=R, 250, w)
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
legend(-15.6,2.4,legend = c("Thinned Estimates", "Standard Normal"), col = c(rgb(116/255,10/255,255/255,alpha=0.7),rgb(43/255,206/255,72/255,alpha=0.7)), lwd = 5, horiz = TRUE, cex = 1.3, seg.len=1, bty = 'n',xpd=NA)
par(op)

