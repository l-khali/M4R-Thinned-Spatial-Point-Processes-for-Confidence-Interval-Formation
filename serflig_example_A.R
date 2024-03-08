library(pracma)
library(latex2exp)

# POISSON PROCESS NOT DISTRIBUTION

poisprocess <- function(lambda, T) {
  arr <- c()
  t <- 0
  while (t < T) {
    inter <- rexp(1, lambda)
    t <- t + inter
    arr <- c(arr, t)
  }
  # inter = rexp(n, lambda)
  # arr = cumsum(inter)
  return(arr)
}


# # using the theoretical results in serfling
# 
# ns <- 10^seq(1,5,0.1)
# lam <- 1
# ders <- c()
# Rs <- c()
# nRs <- c()
# for (n in ns) {
#   print(n)
#   p <- poisprocess(n,lam)
#   #ders <- c(ders, var(p)+mean(p)^2-lam+lam^2-2*lam*mean(p)+2*lam^2)
#   l <- (100000-n)/100000
#   ders <- c(ders, (mean(p)-1/lam)/l)
#   Rs <- c(Rs, -(mean(p)-1/lam)^2)
#   nRs <- c(nRs, -sqrt(n)*(mean(p)-1/lam)^2)
# }
# 
# plot(ns,ders,xlab="log(n)",ylab=TeX("$d_1Mean$"),main="Gateaux derivative of mean")
# plot(ns[-1][-1],Rs[-1][-1],xlab=TeX("$\\log(n)$"),ylab=TeX("$R_n$"),main=TeX("$n$ vs. $R_n$"))
# plot(ns[-1][-1],nRs[-1][-1],xlab=TeX("$\\log(n)$"),ylab=TeX("$\\sqrt{n}R_n$"),main=TeX("$n$ vs. $\\sqrt{n}R_n$"))


# trying my numerical method to see if it matches up

# approaching in the direction of the mean (fixed T)
T <- 100
lambdas <- seq(1,0.000001,-0.01)
ders_mean <- c()

rate <- 1
m <- rate * T

for (lambda in lambdas) {
  Gp <- poisprocess((1-lambda)*rate+lambda*1, T)
  intensity_est <- length(Gp)/T
  print(intensity_est)
  est <- intensity_est * T
  ders_mean <- c(ders_mean,(est-m)/lambda)
}

plot((ders_mean),main="d(mean), approaching in rate", xlab = "lambda", ylab="gataux derivative of mean")

# approaching in direction of T (fixed intensity)
lambdas <- seq(1,0.000001,-0.001)
true_T <- 100
ders_mean2 <- c()
rate <- 1
m <- rate * true_T

for (lambda in lambdas) {
  Gp <- poisprocess(rate, (1-lambda)*true_T+lambda*80)
  est <- rate*max(Gp)
  ders_mean2 <- c(ders_mean2, (est-m)/lambda)
}

plot(ders_mean2,main="d(mean), approaching in T", xlab = "lambda", ylab="gataux derivative of mean")




# POISSON DISTRIBUTION
# approaching in mean
n <- 1000
lambdas <- seq(1,0.000001,-0.01)
ders_mean3 <- c()

mean <- 10

for (lambda in lambdas) {
  Gp <- rpois(n, (1-lambda)*mean+lambda*9)
  ders_mean3 <- c(ders_mean3,(mean(Gp)-mean)/lambda)
}

plot(ders_mean3,main="Limit of gataux derivative numerically, mean", xlab = "lambda", ylab="gataux derivative of mean")


# approaching in n
n <- 10000
lambdas <- seq(1,0.000001,-0.01)
ders_mean4 <- c()

mean <- 10

for (lambda in lambdas) {
  Gp <- rpois((1-lambda)*n + lambda*10, mean)
  ders_mean4 <- c(ders_mean4,(mean(Gp)-mean)/lambda)
}

plot(ders_mean4,main="Limit of gataux derivative numerically, mean", xlab = "lambda", ylab="gataux derivative of mean")


