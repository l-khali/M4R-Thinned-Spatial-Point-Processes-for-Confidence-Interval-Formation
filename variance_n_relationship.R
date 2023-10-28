library(spatstat)

K_var <- function(radius, npoints, window = owin(c(0,1),c(0,1))) {
  # assuming window is unit square
  omega <- area(window)
  beta <- pi * (radius^2) / omega
  gamma <- perimeter(window) * radius / omega
  return(((2*(omega^2)*beta) * (1+0.305*gamma+beta*(-1+0.0132*npoints*gamma)))/(npoints^2))
}

variances <- c()
for (npoints in 10:280) {
  variances <- c(variances, K_var(0.01,npoints))
}

scaling_constant <- function(thinning_param, npoints, radius=seq(0,0.14,0.01), mean=FALSE) {
  if (mean==TRUE) {
    return (mean((K_var(radius,round(npoints*thinning_param))/K_var(radius,npoints))[-1]))
  } else {
    return(K_var(radius,round(npoints*thinning_param))/K_var(radius,npoints))
  }
}

plot(seq(0,0.14,0.01), scaling_constant(0.9, 250), type="l", lty=1, col=1, ylim=c(0,100), xlab="radius", ylab="scaling factor", main="radius vs. scaling factor")
lines(seq(0,0.14,0.01), scaling_constant(0.8, 250), lty=1, col=2)
lines(seq(0,0.14,0.01), scaling_constant(0.7, 250), lty=1, col=3)
lines(seq(0,0.14,0.01), scaling_constant(0.6, 250), lty=1, col=4)
lines(seq(0,0.14,0.01), scaling_constant(0.5, 250), lty=1, col=5)
lines(seq(0,0.14,0.01), scaling_constant(0.4, 250), lty=1, col=6)
lines(seq(0,0.14,0.01), scaling_constant(0.3, 250), lty=1, col=7)
lines(seq(0,0.14,0.01), scaling_constant(0.2, 250), lty=1, col=8)
lines(seq(0,0.14,0.01), scaling_constant(0.1, 250), lty=1, col=9)
legend(0, 95, legend=c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(9,8,7,6,5,4,3,2,1), lty =1, cex=0.7)


plot(seq(0,0.14,0.01), scaling_constant(0.9, 250), type="l", lty=1, col=1, ylim=c(0,5), xlab="radius", ylab="scaling factor", main="radius vs. scaling factor (zoomed in)")
lines(seq(0,0.14,0.01), scaling_constant(0.8, 250), lty=1, col=2)
lines(seq(0,0.14,0.01), scaling_constant(0.7, 250), lty=1, col=3)
lines(seq(0,0.14,0.01), scaling_constant(0.6, 250), lty=1, col=4)
lines(seq(0,0.14,0.01), scaling_constant(0.5, 250), lty=1, col=5)
legend(0, 5, legend=c("0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(5,4,3,2,1), lty =1, cex=0.7)

