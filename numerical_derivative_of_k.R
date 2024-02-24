library(spatstat)
library(pracma)

# attempt 1: moving 1 point in the x-direction, towards a predetermined point 0.5
# using this to approximate derivative and checking limit
# I don't think this is aligned with the actual definition of the Frecht derivative

# will only extract K value for r=0.15
rs <- c(0,0.05,0.1,0.15)
errors <- c()

p <- rpoispp(50, win=owin(c(0,1),c(0,1)))
p2 <- ppp(x=c(p$x, 0.5), y=c(p$y, 0.5), win=owin(c(0,1),c(0,1)))
Kp2 <- Kest(p2,r=rs,correction="isotropic")$iso[4]

new_y <- 0.5
hs <- c(0.01,0.001,0.0001,0.0001,0.00001,0.000001)

for (h in hs) {
  
  Ks <- c()
  xplush <- x + h
  q <- ppp(x=c(p$x, xplush), y=c(p$y, new_y), win=owin(c(0,1),c(0,1)))
  K <- Kest(q,r=rs,correction="isotropic")$iso[4]
  
  new_x <- seq(0.499,0.501,h)
  
  for (x in new_x) {
    q <- ppp(x=c(p$x, x), y=c(p$y, new_y), win=owin(c(0,1),c(0,1)))
    k <- Kest(q,r=rs,correction="isotropic")$iso[4]
    Ks <- c(Ks, k)
  }
  
  poly <- polyfit(new_x, Ks, 1)
  error <- norm(K - Kp2 - poly[1]*xplush - poly[2], type="2")
  errors <- c(errors, error)
}

plot(errors, type="l", main="Numerator")
plot(hs, type="l", main="Denominator")
plot(errors/hs, type="l", main="Ratio")
