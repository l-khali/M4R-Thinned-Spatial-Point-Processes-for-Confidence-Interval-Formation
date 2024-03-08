library(spatstat)

thomas_thinning_simulation <- function(nsim, thinning_param, alpha, kappa, scale, mu) {
  
  # specifying window sizes over which to simulate
  Ws <- seq(1,20,1)
  cover <- rep(c(0),each=length(Ws))
  coverage <- cbind(Ws, cover)
  # indexing radii to select radius of 0.1
  j <- 2
  # actual Ripley's K for Thomas process
  K_actual <- pi*0.1*0.1 + (1/kappa)*(1-exp((-0.1^2)/(4*scale^2)))
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rThomas(kappa, scale, mu, win=square(Ws[i]))
      confidences <- thinning(p, thinning_param, alpha, W=Ws[i],r=c(0.0,0.1))


      if (confidences[j,1] <= K_actual & K_actual <= confidences[j,2]) {
          coverage[i,2] <- coverage[i,2] + 1/nsim
      }
    }
  }
  return(coverage)
}
