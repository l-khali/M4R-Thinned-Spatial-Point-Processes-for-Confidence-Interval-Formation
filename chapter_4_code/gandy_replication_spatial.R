library(pracma)
library(latex2exp)

gandy_asymptotics_spatial <- function(thinning_param, max_W, r = c(0.0,0.1)) {
  # will take n to be in intensity
  Ws <- seq(1,30,1)
  intensity <- 250
  true_k <- pi*r[2]*r[2]
  # intensity <- 10
  covers_thinning <- rep(0, length(Ws))
  for (W in 1:length(Ws)) {
    for (sim in 1:1000){
      print(paste("Current W: ",Ws[W]," Current simulation: ",sim))
      p <- rpoispp(intensity, win=square(W))
      estimates <- c()
      for (bootstrap in 1:100) {
        unif <- runif(npoints(p), 0, 1)
        thinned_p_df <- as.data.frame(p)[which(unif < thinning_param),]
        thinned_p <- as.ppp(thinned_p_df, W=square(W))
        k <- Kest(thinned_p, r=r, correction=c("isotropic"))
        estimates <- c(estimates, k$iso[2])
      }
      conf_interval <- unname(quantile(estimates, probs=c(0.025,0.975), na.rm=TRUE))
      print(conf_interval)
      if (conf_interval[1] <= true_k & true_k <= conf_interval[2]) {
        covers_thinning[W] <- covers_thinning[W] + 1/1000
      }
    }
  } 
  # plot(Ws, covers_thinning)
  # abline(h=0.95,col=2,lty=2)
  return(covers_thinning)
}

k_cover25 <- gandy_asymptotics_spatial(0.25,10)
k_cover5 <- gandy_asymptotics_spatial(0.5,12)
k_cover7 <- gandy_asymptotics_spatial(0.7,10)
k_cover75 <- gandy_asymptotics_spatial(0.75,10)

plot(k_cover25, ylim=c(0,1), type="l")
plot(k_cover5, ylim=c(0,1), type="l")
plot(k_cover7/5, ylim=c(0.5,1), type="l", xlab=TeX("Window Side Length, $W$"), ylab="Cover",main=TeX("Cover of Ripley's $K$ on Homogenous Poisson, $p=0.7$"), cex=1.5)
plot(k_cover75, ylim=c(0,1), type="l")
abline(h=0.95,col=2,lty=2)
lines(k_cover5)
lines(k_cover75)
abline(h=0.95,col=2,lty=2)




