example_140_thinned <- function(thinning_param, nsample = 50, B = 1000, alpha = 0.05) {
  # generating x uniformly on [0,2]
  x <- runif(nsample,0,2)
  x_squared <- x^2
  # generating errors 
  eps <- rnorm(nsample,0,0.2^2)
  # calculating y
  beta0 <- -1
  beta1 <- 2
  beta2 <- -1
  y <- beta0 + beta1 * x + beta2 * x_squared + eps
  
  beta_estimate <- lm(y ~ x + x_squared)
  max_estimate <- - 0.5 * beta_estimate$coefficients[2] / beta_estimate$coefficients[3]
  
  # obtaining estimates of theta by thinned samples
  theta_hat <- c()
  
  for (sim in 1:B) {
    sample <- runif(nsample,0,1)
    thinned_x <- x[which(sample < thinning_param)]
    thinned_x_squared <- thinned_x^2
    thinned_y <- y[which(sample < thinning_param)]
    beta_estimate <- lm(thinned_y ~ thinned_x + thinned_x_squared)
    theta_hat <- c(theta_hat, - 0.5 * beta_estimate$coefficients[2] / beta_estimate$coefficients[3])
  }
  
  sorted_thata_hat <- sqrt(nsample) * (sort(as.numeric(theta_hat)) - max_estimate)
  lower <- max_estimate - sorted_thata_hat[round((B+1)*(1-alpha/2), digits=0)]/sqrt(nsample)
  upper <- max_estimate - sorted_thata_hat[round((B+1)*(alpha/2), digits=0)]/sqrt(nsample)
  
  return(c(lower,upper))
}

example_140_bootstrap <- function(nsample = 50, B = 1000, alpha = 0.05) {
  # generating x uniformly on [0,2]
  x <- runif(nsample,0,2)
  x_squared <- x^2
  # generating errors 
  eps <- rnorm(nsample,0,0.2^2)
  # calculating y
  beta0 <- -1
  beta1 <- 2
  beta2 <- -1
  y <- beta0 + beta1 * x + beta2 * x_squared + eps
  
  beta_estimate <- lm(y ~ x + x_squared)
  max_estimate <- - 0.5 * beta_estimate$coefficients[2] / beta_estimate$coefficients[3]
  
  # obtaining estimates of theta by thinned samples
  theta_hat <- c()
  
  for (sim in 1:B) {
    sample <- sample(1:nsample, 50, replace=TRUE)
    x_star <- x[sample]
    x_star_squared <- x_star^2
    y_star <- y[sample]
    beta_estimate <- lm(y_star ~ x_star + x_star_squared)
    theta_hat <- c(theta_hat, - 0.5 * beta_estimate$coefficients[2] / beta_estimate$coefficients[3])
  }
  
  sorted_thata_hat <- sqrt(nsample) * (sort(as.numeric(theta_hat)) - max_estimate)
  lower <- max_estimate - sorted_thata_hat[round((B+1)*(1-alpha/2), digits=0)]/sqrt(nsample)
  upper <- max_estimate - sorted_thata_hat[round((B+1)*(alpha/2), digits=0)]/sqrt(nsample)
  
  return(c(lower,upper))
}

coverage_example_140_thinning <- function(nsim=1000, alpha=0.05) {
  thinning_params <- seq(0.3,1,0.05)
  coverages <- c()
  for (thinning_param in thinning_params) {
    coverage <- 0
    for (sim in 1:nsim) {
      interval <- example_140_thinned(thinning_param, alpha=alpha)
      # 1 is the true value of theta
      if (interval[1] <= 1 & 1<=interval[2]) {
        coverage = coverage + 1/nsim
      }
    }
    coverages <- c(coverages, coverage)
  }
  return(coverages)
}

coverage_example_140_bootstrap <- function(nsim=1000, alpha=0.05) {
  coverage <- 0
  for (sim in 1:nsim) {
    interval <- example_140_bootstrap(alpha=alpha)
    print(interval)
      # 1 is the true value of theta
    if (interval[1] <= 1 & 1<=interval[2]) {
      coverage = coverage + 1/nsim
    }
  }
  return(coverage)
}

example_140_coverage <- coverage_example_140_thinning()
plot(seq(0.3,1,0.05), example_140_coverage, type="l", main="Example 104 thinned")
abline(h = 0.95, col=2, lty=2)

example_140_bs_coverage <- coverage_example_140_bootstrap()
print(paste0("Coverage using boostrap method = ", example_140_bs_coverage))
