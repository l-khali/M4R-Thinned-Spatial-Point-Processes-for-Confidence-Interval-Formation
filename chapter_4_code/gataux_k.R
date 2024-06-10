library(spatstat)

# For a fixed radius, the value of ripley's K is poisson distributed with 
# mean (lambda*pi*r^2).
# Attempt to find gatuax derivative by showing mean of poisson is gataux

# approaching in lambda
mus <- 10^(-seq(-1,8,0.01))
lambda <- 10
lambda_G <- 9
ders <- c()
diffs <- c()
true <- lambda_G-lambda

for (mu in mus) {
  Gp <- rpois(5000, ((1-mu)*lambda + mu*lambda_G))
  diffs <- c(diffs, (mean(Gp)-lambda))
  ders <- c(ders, (mean(Gp)-lambda)/mu)
}

plot(ders)
plot(diffs)


# approaching in n
mus <- seq(1,0.000001,-0.01)
lambda <- 10
ders <- c()
diffs <- c()

for (mu in mus) {
  Gp <- rpois((1-mu)*1000 + mu*10, lambda)
  ders <- c(ders, (mean(Gp)-lambda)/mu)
  diffs <- c(diffs, mean(Gp)-lambda)
}

plot(ders)
plot(diffs)


# variance: approaching in lambda
mus <- seq(1,0.000001,-0.01)
lambda <- 10
ders <- c()
diffs <- c()
true <- 9+9^2-10+10^2-2*9*10+2*10^2

for (mu in mus) {
  Gp <- rpois(5000, ((1-mu)*lambda + mu*9))
  diffs <- c(diffs, (var(Gp)-lambda))
  ders <- c(ders, (var(Gp)-lambda)/mu)
}

plot(ders)
