# do this with ripleys k for now
# need to adapt the thinning function to return the normalised statistics

thinning_jh <- function(data, thinning_param, alpha, R=99, radius = c(0,0.1)) {
  # confidence intervals calculated using quantiles of samples
  process_df <- as.data.frame(data)
  npoints <- nrow(process_df)
  k_vals <- c()
  K_est <- as.data.frame(Kest(data, r = radius, correction=c("isotropic")))["iso"]
  K_est <- K_est["iso"][2,1]
  js <- c()
  for (i in 1:R) {
    # unif <- runif(npoints, 0, 1)
    # subprocess_df <- process_df[which(unif < thinning_param),]
    # subprocess <- as.ppp(subprocess_df, owin(c(0,1),c(0,1)))
    subprocess <- rthin(data, thinning_param)
    rs <- npoints(subprocess)
    ds <- npoints - rs
    k <- Kest(subprocess, r=radius, correction=c("isotropic"))
    k_val <- as.data.frame(k)["iso"]
    k_vals <- c(k_vals, k_val["iso"][2,1])
    js <- c(js, (k_val["iso"][2,1] - K_est)*sqrt(npoints * rs/ds))
  }
  k_vals[is.nan(k_vals)] <- 0 # subprocesses can have 0 points with low thinning param
  # js <- js / var(as.numeric(k_vals))
  return(js)
}

# grid of plots
par(mfrow = c(3, 3))
op <- par(mfrow = c(3,3),
          oma = c(5,4,2,1) + 0.1,
          mar = c(2,2,2,2) + 0.1)

# Homogenous Poisson Process (independent)
for (thinning_param in c(0.25,0.5,0.75)) {
  R <- 2000
  ns <- c(2.3,4,150)
  for (n in ns) {
    # homogenous poisson with mean 50
    p <- rpoispp(5*n)
    npoint <- nrow(as.data.frame(p))
    # obtaining jackknife histogram
    thinned_estimates <- thinning_jh(p, thinning_param, 0.05, R=R)
    normal_dist <- rnorm(R, 0, 1)
    print(paste("n", n,"length",length(p), "estimates",length(thinned_estimates)))
    # qqplot to observe convergence in distribution
    qqnorm(thinned_estimates, main=paste("p =", thinning_param, ", n =", npoint), ylab = "", xlab="")
    qqline(thinned_estimates)
    # qqplot(normal_dist, thinned_estimates, main=paste("p =", thinning_param, ", n =", npoint), ylab = "", xlab="")
    # qqline(thinned_estimates, distribution= function(p) qnorm(p, 0, 1), col = 2)
  }
}
title("Homogenous Poisson Point Process",line=0.5,outer=TRUE, cex.main=1.5)
title(xlab="Standard Normal", ylab="Normalised Thinned Estimates", outer=TRUE, line=1, cex.lab=1.5)
par(op)


# new grid of plots
par(mfrow = c(3, 3))
op <- par(mfrow = c(3,3),
          oma = c(5,4,2,1) + 0.1,
          mar = c(2,2,2,2) + 0.1)

# Neyman-Scott Process (independent)
for (thinning_param in c(0.25,0.5,0.75)) {
  R <- 2000
  ns <- c(1.7,3,15)
  for (n in ns) {
    # neymann scott
    set.seed(2)
    p <- rMatClust(2.5*n,0.1,1*n, win=owin(c(0,1),c(0,1)))
    npoint <- nrow(as.data.frame(p))
    print(paste("number of points", npoint))
    # obtaining jackknife histogram
    thinned_estimates <- thinning_jh(p, thinning_param, 0.05, R=R)
    normal_dist <- rnorm(R, 0, 1)
    # qqplot to observe convergence in distribution
    qqplot(normal_dist, thinned_estimates, main=paste("p =", thinning_param, ", n =", npoint), ylab = "", xlab="")
    qqline(thinned_estimates, distribution= function(p) qnorm(p, 0, 1), col = 2)
  }
}
title("Neyman-Scott Point Process",line=1,outer=TRUE, cex.main=1.5)
title(xlab="Standard Normal", ylab="Normalised Thinned Estimates", outer=TRUE, line=1, cex.lab=1.5)
par(op)

