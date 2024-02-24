# take true distribution, F, to be hom poisson with n=10000 
# take G to be empirical distribution varying over n=9900,...,10000

rs <- c(0,0.05,0.1,0.15)
lambdas <- seq(1,0,-0.001)
ders <- c()

Fp <- rpoispp(250)
true_K <- Kest(Fp,r=rs,correction="isotropic")$iso[4]

for (lambda in lambdas) {
  print(250+lambda*(249.999-250))
  Gp <- rpoispp(250+lambda*(249.999-250))
  ders <- c(ders,(Kest(Gp,r=rs,correction="isotropic")$iso[4]-true_K)/lambda)
}

plot(ders[-100],main="Limit of gataux derivative (seems to not exist)", xlab = "lambda", ylab="gataux derivative of K", type="l")

# testing if sample variance is gataux differentiable in the same way
# this shouldnt diverge since variance is differentiable (Serfling)

ders_mean <- c()

Fp <- rpois(100, 250)
true_mean <- mean(Fp)

for (lambda in lambdas) {
  print(250+lambda*(249-250))
  Gp <- rpois(250+lambda*(249-250), 250)
  m <- mean(rpois(250)) + lambda*mean(rpois(249))
  ders_mean <- c(ders_mean,(mean(Gp)-true_mean)/lambda)
}

plot((ders_mean[-100]),main="Limit of gataux derivative", xlab = "lambda", ylab="gataux derivative of mean")
abline(h=mean(Gp)-mean(Fp),col=2)
