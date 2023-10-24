library(spatstat)

thinning <- function(data, thinning_param, alpha, R=99) {
  process_df <- as.data.frame(data)
  npoints <- nrow(process_df)
  r <- seq(0.0, 0.14, 0.01)
  k_vals <- data.frame(r)
  for (i in 1:R) {
    unif <- runif(npoints, 0, 1)
    subprocess_df <- process_df[which(unif < thinning_param),]
    subprocess <- as.ppp(subprocess_df, owin(c(0,1),c(0,1)))
    k <- Kest(subprocess, r=r, correction=c("isotropic"))
    k_vals <- cbind(k_vals, as.data.frame(k)["iso"])
  }
  K_est <- as.data.frame(Kest(data, r=r, correction=c("isotropic")))["iso"]
  lower_approx <- c()
  upper_approx <- c()
  
  for (radius in 1:15) {
    sorted_k_vals <- sort(as.numeric(k_vals[radius,-1]))
    lower <- sorted_k_vals[(R+1)*(1-alpha/2)]
    upper <- sorted_k_vals[(R+1)*(alpha/2)]
    lower_approx <- c(lower_approx, 2*K_est[radius,] - lower)
    upper_approx <- c(upper_approx, 2*K_est[radius,] - upper)
  }
  return(cbind(lower_approx, upper_approx)) 
}

thinning_1 <- poisson_simulation_thinning(500, 250, 1, 0.05)
plot(thinning_1[-1,], type="l", lty=1, ylim=c(0,1), main="Poisson: thinning (1.0)")
thinning_09 <- poisson_simulation_thinning(500, 250, 0.9, 0.05)
plot(thinning_09[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.9)")
thinning_08 <- poisson_simulation_thinning(500, 250, 0.8, 0.05)
plot(thinning_08[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.8)")
thinning_07 <- poisson_simulation_thinning(500, 250, 0.7, 0.05)
plot(thinning_07[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.7)")
thinning_06 <- poisson_simulation_thinning(500, 250, 0.6, 0.05)
plot(thinning_06[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.6)")
thinning_05 <- poisson_simulation_thinning(500, 250, 0.5, 0.05)
plot(thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.5)")
thinning_04 <- poisson_simulation_thinning(500, 250, 0.4, 0.05)
plot(thinning_04[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.4)")
thinning_03 <- poisson_simulation_thinning(500, 250, 0.3, 0.05)
plot(thinning_03[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.3)")
thinning_02 <- poisson_simulation_thinning(500, 250, 0.2, 0.05)
plot(thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.2)")
thinning_01 <- poisson_simulation_thinning(500, 250, 0.1, 0.05)
plot(thinning_01[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning (0.1)")

plot(thinning_09[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=1)
lines(thinning_08[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=2)
lines(thinning_07[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=3)
lines(thinning_06[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=4)
lines(thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=5)
lines(thinning_04[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=6)
lines(thinning_03[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=7)
lines(thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=8)
lines(thinning_01[-1,], type="l", lty=1, ylim=c(0.5,1), main="Poisson: thinning", col=9)
legend(0.12, 0.75, legend=c("0.9", "0.8", "0.7", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1"),
       col=c(1,2,3,4,5,6,7,8,9), lty =1, cex=0.8)
