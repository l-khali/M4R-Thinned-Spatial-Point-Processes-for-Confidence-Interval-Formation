library("GiniDistance")
library("hawkesbow")

# FIND A STATISTIC WITH A HIGHER ORDER REMAINEDER TERM


# Observing how well thinning method works asymptotically on one dimensional processes
# estimate variance
ests <- c()
for (sim in 1:10000) {
  p <- rpois(1000,50)
  ests <- c(ests, gmd(p))
}
v <- var(ests)


# Thinning function for one dimensional process which statistic being Gini mean difference
thinning1d_gmd_jh <- function(data, thinning_param, alpha, R=99) {
  n <- length(data)
  process_df <- as.data.frame(data)
  npoints <- nrow(process_df)
  gmds <- c()
  js <- c()
  for (i in 1:R) {
    # thinning the given process
    unif <- runif(npoints, 0, 1)
    subprocess_df <- process_df[which(unif < thinning_param),]
    r <- length(as.numeric(subprocess_df))
    d <- npoints - r
    # calculating jackknife histogram
    j <- (gmd(as.numeric(subprocess_df)) - gmd(as.numeric(data)))*sqrt(n * r/d)
    gmds <- c(gmds, gmd(as.numeric(subprocess_df)))
    js <- c(js, j)
  }
  gmds[is.nan(gmds)] <- 0 # subprocesses can have 0 points with low thinning param
  print(var(gmds))
  js <- js / v
  return(js)
}

# grid of plots
par(mfrow = c(3, 3))
op <- par(mfrow = c(3,3),
          oma = c(5,4,2,1) + 0.1,
          mar = c(2,2,2,2) + 0.1)

# Homogenous Poisson Process (independent)
for (thinning_param in c(0.2,0.5,0.8)) {
  R <- 2000
  ns <- c(5,20,100)
  for (n in ns) {
    # homogenous poisson with mean 50
    p <- rpois(n,50)
    npoint <- nrow(as.data.frame(p))
    # obtaining jackknife histogram
    thinned_estimates <- thinning1d_gmd_jh(p, thinning_param, 0.05, R=R)
    normal_dist <- rnorm(R, 0, 1)
    # qqplot to observe convergence in distribution
    qqplot(normal_dist, thinned_estimates, main=paste("p =", thinning_param, ", n =", npoint), ylab = "", xlab="")
    qqline(thinned_estimates, distribution= function(p) qnorm(p, 0, 1), col = 2)
  }
}
title("Homogenous Poisson Process",line=0.5,outer=TRUE, cex.main=2)
title(xlab="Standard Normal", ylab="Normalised Thinned Estimates", outer=TRUE, line=1, cex.lab=2)
par(op)

# Correlated Hawkes process
# grid of plots
# grid of plots
par(mfrow = c(3, 3))
op <- par(mfrow = c(3,3),
          oma = c(5,4,2,1) + 0.1,
          mar = c(2,2,2,2) + 0.1)

intensity <- function(x) {
  return(100 - abs(x-100))
}

for (thinning_param in c(0.2,0.5,0.8)) {
  R <- 2000
  ns <- c(5,50,500000)
  for (n in ns) {
    # homogenous poisson with mean 50
    set.seed(1)
    p <- hawkes(10, fun=log10(n), repr=0.1*log10(n), family="exp", rate=2)$p
    # obtaining jackknife histogram
    thinned_estimates <- thinning1d_gmd_jh(p, thinning_param, 0.05, R=R)
    normal_dist <- rnorm(R, 0, 1)
    # qqplot to observe convergence in distribution
    qqplot(normal_dist, thinned_estimates, main=paste("p =", thinning_param, ", n =", length(p)), ylab = "", xlab="")
    qqline(thinned_estimates, distribution= function(p) qnorm(p, 0, 1), col = 2)
  }
}
title("Hawkes Process",line=0.5,outer=TRUE, cex.main=2)
title(xlab="Standard Normal", ylab="Normalised Thinned Estimates", outer=TRUE, line=1, cex.lab=2)
par(op)


