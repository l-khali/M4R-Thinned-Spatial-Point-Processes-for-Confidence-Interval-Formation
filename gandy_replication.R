ns <- c(10,20,40,80,160,320)
intensity <- 10
covers <- rep(0, length(ns))
for (n in 1:length(ns)) {
  for (sim in 1:2000){
    p <- rexp(ns[n], intensity)
    estimates <- c()
    for (bootstrap in 1:500) {
      p_bootstrapped <- sample(p,size=ns[n],replace=TRUE)
      estimates <- c(estimates, mean(p_bootstrapped))
    }
    conf_interval <- unname(quantile(estimates, probs=c(0.025,0.975)))
    print(conf_interval)
    if (conf_interval[1] <= 1/intensity & 1/intensity <= conf_interval[2]) {
      covers[n] <- covers[n] + 1/2000
    }
  }
} 
plot(covers)



gandy_asymptotics <- function(thinning_param, max_n) {
  # ns <- c(10,20,40,80,160,320,640,1280,2560,5120,10240)
  ns <- 10*(2^0:max_n)
  intensity <- 10
  covers_thinning <- rep(0, length(ns))
  for (n in 1:length(ns)) {
    for (sim in 1:2000){
      print(paste("Current n: ",ns[n]," Current simulation: ",sim))
      p <- rexp(ns[n], intensity)
      estimates <- c()
      for (bootstrap in 1:500) {
        unif <- runif(ns[n], 0, 1)
        p_thinned <- p[which(unif < thinning_param)]
        estimates <- c(estimates, mean(p_thinned))
      }
      conf_interval <- unname(quantile(estimates, probs=c(0.025,0.975), na.rm=TRUE))
      if (conf_interval[1] <= 1/intensity & 1/intensity <= conf_interval[2]) {
        covers_thinning[n] <- covers_thinning[n] + 1/2000
      }
    }
  } 
  plot(ns, covers_thinning)
  abline(h=0.95,col=2,lty=2)
  return(covers_thinning)
}

exp_cover25 <- gandy_asymptotics(0.25,10)
exp_cover5 <- gandy_asymptotics(0.5,12)
exp_cover75 <- gandy_asymptotics(0.75,10)

plot(exp_cover25, ylim=c(0,1),type="l")
lines(exp_cover5)
lines(exp_cover75)
abline(h=0.95,col=2,lty=2)





