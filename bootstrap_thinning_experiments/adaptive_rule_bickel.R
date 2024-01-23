adaptive_rule_normal <- function() {
  mjs <- 0.9 ^ (1:20)
  rho <- c()
  for (j in 1:(length(mjs)-1)) {
    print(paste0("Current mj:", mjs[j]))
    temp1 <- sort(normal_thinning(thinning_param = mjs[j], B=10000))
    temp2 <- sort(normal_thinning(thinning_param = mjs[j+1], B=10000))
    rho <- c(rho, max(abs(temp1 - temp2)))
  }
  i <- which.min(rho)
  m_hat <- mjs[i]
  return(m_hat)
}

adaptive_rule_exponential <- function() {
  mjs <- 0.9 ^ (1:20)
  rho <- c()
  for (j in 1:(length(mjs)-1)) {
    print(paste0("Current mj:", mjs[j]))
    temp1 <- sort(exponential_experiment(thinning_param = mjs[j], B=10000))
    temp2 <- sort(exponential_experiment(thinning_param = mjs[j+1], B=10000))
    rho <- c(rho, max(abs(temp1 - temp2)))
  }
  i <- which.min(rho)
  m_hat <- mjs[i]
  return(m_hat)
}

adaptive_rule_poisson <- function() {
  mjs <- 0.9 ^ (1:20)
  rho <- c()
  for (j in 1:(length(mjs)-1)) {
    print(paste0("Current mj:", mjs[j]))
    temp1 <- sort(poisson_experiment(thinning_param = mjs[j], B=10000))
    temp2 <- sort(poisson_experiment(thinning_param = mjs[j+1], B=10000))
    rho <- c(rho, max(abs(temp1 - temp2)))
  }
  i <- which.min(rho)
  m_hat <- mjs[i]
  return(m_hat)
}

adaptive_rule_exp_median <- function() {
  mjs <- 0.9 ^ (1:20)
  rho <- c()
  for (j in 1:(length(mjs)-1)) {
    print(paste0("Current mj:", mjs[j]))
    temp1 <- sort(efron_2(thinning_param = mjs[j], B=10000))
    temp2 <- sort(efron_2(thinning_param = mjs[j+1], B=10000))
    rho <- c(rho, max(abs(temp1 - temp2)))
  }
  i <- which.min(rho)
  m_hat <- mjs[i]
  return(m_hat)
}

normal_ms <- c()
exponential_ms <- c()
poisson_ms <- c()

for (i in 1:50) {
  normal_ms <- c(normal_ms, adaptive_rule_normal())
  exponential_ms <- c(exponential_ms, adaptive_rule_exponential())
  poisson_ms <- c(poisson_ms, adaptive_rule_poisson())
}

print(paste0("Normal:", mean(normal_ms)))
print(paste0("Exponential:", mean(exponential_ms)))
print(paste0("Poisson:", mean(poisson_ms)))

efron_ms <- c()
for (i in 1:50) {
  efron_ms <- c(efron_ms, adaptive_rule_exp_median())
}

print(paste0("Efron 2:", mean(efron_ms)))