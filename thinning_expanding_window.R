library(spatstat)

# using thinning to find confidence intervals in the limit 
# as window size grows (and so n grows)

poisson_expanding_window <- function(nsim, thinning_param, alpha, intensity) {
  
  # specifying window sizes over which to simulate
  Ws <- seq(0.5,2,0.2)
  cover <- rep(c(0),each=length(Ws))
  coverage <- cbind(Ws, cover)
  # indexing radii to select radius of 0.1
  j <- 2
  # actual Ripley's K for Thomas process
  K_actual <- pi*0.1*0.1
  
  for (i in 1:length(Ws)) {
    print(paste0("Current window length:", Ws[i]))
    for (sim in 1:nsim) {
      p <- rpoispp(intensity, win=square(Ws[i]))
      confidences <- thinning(p, thinning_param, alpha, W=Ws[i],r=c(0.0,0.1))
      print(confidences)
      
      if (confidences[j,1] <= K_actual & K_actual <= confidences[j,2]) {
        coverage[i,2] <- coverage[i,2] + 1/nsim
      }
    }
  }
  print(coverage)
  return(coverage)
}

cover1 <- poisson_expanding_window(200,0.1,0.05,250) # some thinned samples end up having no points and cause errors
cover2 <- poisson_expanding_window(200,0.2,0.05,250)
cover3 <- poisson_expanding_window(200,0.3,0.05,250)
cover4 <- poisson_expanding_window(200,0.4,0.05,250)
cover5 <- poisson_expanding_window(200,0.5,0.05,250)
cover6 <- poisson_expanding_window(200,0.6,0.05,250)
cover7 <- poisson_expanding_window(200,0.7,0.05,250)
cover8 <- poisson_expanding_window(200,0.8,0.05,250)
cover9 <- poisson_expanding_window(200,0.9,0.05,250)

plot(cover2)
lines(cover3)
lines(cover4)
lines(cover5)
lines(cover6)
lines(cover7)
lines(cover8)
lines(cover9)

