adaptive_rule_normal <- function() {
  mjs <- 0.9 ^ (1:40)
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