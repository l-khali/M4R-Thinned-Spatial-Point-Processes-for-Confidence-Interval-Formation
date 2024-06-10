ns <- seq(from = 10, to = 4000, by = 200)
rs <- c(0,0.05,0.1,0.15)
Rs <- data.frame(rs)
Vs <- data.frame(rs)
# Ds <- data.frame(rs)
for (n in ns) {
  Ravg <- data.frame(c(0,0,0,0))
  Vavg <- data.frame(c(0,0,0,0))
  Davg <- data.frame(c(0,0,0,0))
  for (i in 1:5) {
    p <- rpoispp(n, win=owin(c(0,1),c(0,1)))
    R <- sqrt(n) * (Kest(p,r=rs,correction="isotropic")$iso - pi*rs^2 - 2*pi*rs*pcf(p,r=rs)$iso + 2*pi*rs)
    V <- sqrt(n) * (Kest(p,r=rs,correction="isotropic")$iso - pi*rs^2)
    # D <- sqrt(n) * Kest(p,r=rs,correction="isotropic")$iso
    # k <- Kest(p,r=rs,correction="isotropic")$iso
    # R <- sqrt(n) * (k - pi*rs^2 - k*(2*pi*rs*pcf(p,r=rs)$iso)+ 2*pi*rs)
    Ravg <- Ravg + R/5
    Vavg <- Vavg + V/5
    # Davg <- Davg + D/5
  }
  Rs <- cbind(Rs,Ravg)
  Vs <- cbind(Vs,Vavg)
  # Ds <- cbind(Ds,Davg)
}

plot(as.numeric(Rs[2,-1]))
abline(h=0,col=2)
plot(as.numeric(Rs[3,-1]),col=2)
abline(h=0,col=2)
plot(as.numeric(Rs[4,-1]),col=3, xlab="seq(from = 10, to = 20000, by = 500)")
abline(h=0,col=2)

# plot(as.numeric(Vs[2,-1]))
# abline(h=0,col=2)
# plot(as.numeric(Vs[3,-1]),col=2)
# abline(h=0,col=2)
# plot(as.numeric(Vs[4,-1]),col=3, xlab="seq(from = 10, to = 40000, by = 2000)")
# abline(h=0,col=2)

# plot(as.numeric(Ds[2,-1]))
# abline(h=0,col=2)
# plot(as.numeric(Ds[3,-1]),col=2)
# abline(h=0,col=2)
# plot(as.numeric(Ds[4,-1]),col=3, xlab="seq(from = 10, to = 20000, by = 500)")
# abline(h=0,col=2)


                                        