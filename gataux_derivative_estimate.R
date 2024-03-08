library(spatstat)
# take true distribution, F, to be hom poisson with n=10000 
# take G to be empirical distribution varying over n=9900,...,10000

# ripleys K
rs <- c(0,0.01,0.1)
lambdas <- seq(1,0.000001,-0.001)
ders <- c()

# Fp <- rpoispp(250)
# true_K <- Kest(Fp,r=rs,correction="isotropic")$iso[4]
true_K <- pi*0.1*0.1

for (lambda in lambdas) {
  Gp <- rpoispp((1-lambda)*10000+lambda*(100))
  k <- Kest(Gp,r=rs,correction="isotropic")$iso[3]
  ders <- c(ders,(k-true_K)/lambda)
}

plot(ders,main="Limit of gataux derivative of K, lambda", xlab = "lambda", ylab="gataux derivative of K")


# gradually making region of process larger instead
rs <- c(0,0.01,0.1)
lambdas <- seq(1,0.000001,-0.01)
ders <- c()
diffs <- c()

# Fp <- rpoispp(250)
# true_K <- Kest(Fp,r=rs,correction="isotropic")$iso[4]
true_K <- pi*0.1*0.1

for (lambda in lambdas) {
  Gp <- rpoispp(5, win=owin(c(0,(1-lambda)*100+lambda*1),c(0,(1-lambda)*100+lambda*1)))
  print(Gp)
  k <- Kest(Gp,r=rs,correction="isotropic")$iso[3]
  diffs <- c(diffs, k-true_K)
  ders <- c(ders,(k-true_K)/lambda)
}

plot(ders,main="Limit of gataux derivative of K, window", xlab = "lambda", ylab="gataux derivative of K")  
plot(diffs)


# # testing if sample variance is gataux differentiable in the same way
# # this shouldnt diverge since variance is differentiable (Serfling)
# 
# ders_mean <- c()
# 
# Fp <- rpois(100, 250)
# true_mean <- mean(Fp)
# 
# for (lambda in lambdas) {
#   Gp <- rpois(100, 250+lambda*(249-250))
#   ders_mean <- c(ders_mean,(mean(Gp)-true_mean)/lambda)
# }
# 
# plot((ders_mean[-100]),main="Limit of gataux derivative", xlab = "lambda", ylab="gataux derivative of mean")
# abline(h=mean(rpois(100,249))-mean(Fp),col=2)


# intensity estimate instead of ripleys
rs <- c(0,0.01,0.1)
lambdas <- seq(1,0.000001,-0.01)
ders <- c()
diffs <- c()

# Fp <- rpoispp(250)
# true_K <- Kest(Fp,r=rs,correction="isotropic")$iso[4]
true_intensity <- 50

for (lambda in lambdas) {
  Gp <- rpoispp(50, win=owin(c(0,(1-lambda)*10+lambda*0.1),c(0,(1-lambda)*10+lambda*0.1)))
  i <- npoints(Gp)/((1-lambda)*10+lambda*0.1)^2
  diffs <- c(diffs, i-true_intensity)
  ders <- c(ders,(i-true_intensity)/lambda)
}

plot(ders,main="Limit of gataux derivative (seems to not exist)", xlab = "lambda", ylab="gataux derivative of intensity",type="l")  

