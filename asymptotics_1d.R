library("GiniDistance")
source("thinning.R")

# Observing how well thinning method works asymptotically on one dimensional processes

# Thinning function for one dimensional process which statistic being mean (DO GINI MEAN DIFFERENCE)
thinning1d_mean <- function(data, thinning_param, alpha, R=99) {
  # confidence intervals calculated using quantiles of samples
  process_df <- as.data.frame(data)
  npoints <- nrow(process_df)
  gmds <- c()
  for (i in 1:R) {
    unif <- runif(npoints, 0, 1)
    subprocess_df <- process_df[which(unif < thinning_param),]
    gmds <- cbind(gmds, mean(as.numeric(subprocess_df)))
  }
  gmd_est <- mean(as.numeric(data))
  
  sorted_vals <- sort(as.numeric(gmds))
  lower <- sorted_vals[(R+1)*(1-alpha/2)]
  upper <- sorted_vals[(R+1)*(alpha/2)]
  lower_approx <- 2*gmd_est - lower
  upper_approx <- 2*gmd_est - upper
  
  return(cbind(lower_approx, upper_approx)) 
}

# Homogenous Poisson Process (independent)
coverage <- c()
for (n in c(10,100,1000,10000)) {
  nsim <- 100
  for (i in 1:nsim) {
    print(paste0("Current simulation:",i))
    
    # generating homogenous poisson process
    p <- rpois(n, 100)
    cover <- 0
    
    conf_interval <- thinning1d_mean(p, 0.5, 0.05)
    if (conf_interval[1] <= 100 & 100 <= conf_interval[2]) {
      cover <- cover + 1/nsim
      print("covered!")
    }
  }
  coverage <- cbind(coverage, cover)
}

plot(coverage)

# Inhomogenous Poisson Process (independent conditional on intensity)

# some very correlated process

