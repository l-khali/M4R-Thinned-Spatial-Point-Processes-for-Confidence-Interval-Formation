library(spatstat)
library(pracma)
library(latex2exp)

# attempt 2: using limit as n goes to infinity so that F_n goes to F

# maxns <- c(10,100,500,1000,5000,10000,50000,100000, 500000)
maxns <- seq(9900,10000)
errors <- c()
second_ders <- c()

for (maxn in maxns) {
  rs <- c(0,0.05,0.1,0.15)
  ns <- seq(9*maxn/10,11*maxn/10,maxn/100)
  Ks <- c()
  
  for (n in ns) {
    p <- rpoispp(n, win=owin(c(0,1),c(0,1)))
    Ks <- c(Ks, Kest(p,r=rs,correction="isotropic")$iso[4])
  }
  
  poly <- polyfit(ns, Ks, 2)
  
  Kest <- Kest(p,r=rs,correction="isotropic")$iso[4]
  #error <- norm(Kest - poly[1]*ns - poly[2], type="2")
  error <- norm(Kest - poly[1]*ns - pi*0.15*0.15, type="2")
  errors <- c(errors, error)
  second_der <- poly[2]*ns
  second_ders <- c(second_ders, norm(second_der))
}

# plot(maxns[-1][-1], errors[-1][-1], main=TeX("$n$ vs. $R_n$"), xlab=TeX("$n$"), ylab=TeX("$R_n$"))
# plot(maxns[-1][-1], (sqrt(maxns) * errors)[-1][-1], main=TeX("$n$ vs. $\\sqrt{n}R_n$"), xlab=TeX("$n$"), ylab=TeX("$\\sqrt{n}R_n$"))
plot(maxns, (sqrt(maxns) * second_ders), main=TeX("$n$ vs. $\\sqrt{n}R_n$"), xlab=TeX("$n$"), ylab=TeX("$\\sqrt{n}R_n$"))



