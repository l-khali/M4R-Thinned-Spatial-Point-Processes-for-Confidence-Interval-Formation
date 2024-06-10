library(mvtnorm)
library(MASS)

efron_3 <- function(thinning_param = 0.5, B = 1000, n = 50, m = 50) {
  X_population <- rmvnorm(m, mean = c(-0.5, 0), sigma = matrix(c(1, 0, 0, 1), nrow = 2))
  Y_population <- rmvnorm(n, mean = c(0.5, 0), sigma = matrix(c(1, 0, 0, 1), nrow = 2))
  
  bootstrapped_Rs <- c()
  thinned_Rs <- c()
  
  for (b in 1:B) {
    bootstrap_sample_x <- X_population[sample(nrow(X_population), m), ]
    bootstrap_sample_y <- Y_population[sample(nrow(Y_population), n), ]
    x_star_bar_boot <- matrix(c(mean(bootstrap_sample_x[,1]),mean(bootstrap_sample_x[,2])))
    y_star_bar_boot <- matrix(c(mean(bootstrap_sample_y[,1]),mean(bootstrap_sample_y[,2])))
    
    rand <- runif(m,0,1)
    thinned_sample_x <- X_population[which(rand < thinning_param),]
    rand <- runif(n,0,1)
    thinned_sample_y <- Y_population[which(rand < thinning_param),]
    x_star_bar_thin <- matrix(c(mean(thinned_sample_x[,1]),mean(thinned_sample_x[,2])))
    y_star_bar_thin <- matrix(c(mean(thinned_sample_y[,1]),mean(thinned_sample_y[,2])))
    
    diff_x_boot <- matrix(bootstrap_sample_x, ncol=2) - matrix(x_star_bar_boot,ncol=2,nrow=m,byrow=TRUE)
    diff_y_boot <- matrix(bootstrap_sample_y, ncol=2) - matrix(y_star_bar_boot,ncol=2,nrow=n,byrow=TRUE)
    S_boot <- sum(diag(t(diff_x_boot) %*% diff_x_boot)) + sum(diag(t(diff_y_boot) %*% diff_y_boot))
    
    diff_x_thin <- matrix(bootstrap_sample_x, ncol=2) - matrix(x_star_bar_boot,ncol=2,nrow=m,byrow=TRUE)
    diff_y_thin <- matrix(bootstrap_sample_y, ncol=2) - matrix(y_star_bar_boot,ncol=2,nrow=n,byrow=TRUE)
    S_thin <- sum(diag(t(diff_x_thin) %*% diff_x_thin)) + sum(diag(t(diff_y_thin) %*% diff_y_thin))
    
    error_star_boot <- 0
    error_star_thin <- 0
    error_boot <- 0
    error_thin <- 0
    
    for (x_star in 1:nrow(thinned_sample_x)) {
      xi <- matrix(thinned_sample_x[x_star,])
      if ((t(y_star_bar_thin - x_star_bar_thin) %*% (xi - (x_star_bar_thin + y_star_bar_thin)/2)) > log(m/n)) { 
        error_star_thin <- error_star_thin + 1/m}
      if ((t(y_star_bar_boot - x_star_bar_boot) %*% (xi - (x_star_bar_boot + y_star_bar_boot)/2)) > log(m/n)) { 
        error_star_boot <- error_star_boot + 1/m}
    }
    
    for (x_star in 1:nrow(X_population)) {
      xi <- matrix(X_population[x_star,])
      if ((t(y_star_bar_thin - x_star_bar_thin) %*% (xi - (x_star_bar_thin + y_star_bar_thin)/2)) > log(m/n)) { 
        error_thin <- error_thin + 1/m}
      if ((t(y_star_bar_boot - x_star_bar_boot) %*% (xi - (x_star_bar_boot + y_star_bar_boot)/2)) > log(m/n)) { 
        error_boot <- error_boot + 1/m}
    }
    
    bootstrapped_Rs <- c(bootstrapped_Rs, error_boot - error_star_boot)
    thinned_Rs <- c(thinned_Rs, error_thin - error_star_thin)
    
  }
  qqplot(bootstrapped_Rs, thinned_Rs, main=paste("Error of discriminant analysis, thinning parameter =", thinning_param))
  abline(0, 1, col = 'red')
}

par(mfrow = c(2, 5))
for (param in seq(0.2,1,0.1)) {
  efron_3(thinning_param = param, B=1000, n=50)
}
