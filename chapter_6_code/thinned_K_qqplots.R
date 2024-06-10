library(spatstat)
library(latex2exp)


k_thinned_estimates <- function(thinning_param) {
  thinned_k_hats <- c()
  normalised_thinned_k_hats <- c()
  for (i in 1:1) {
    p <- rpoispp(10000,win=square(1))
    k_hat <- Kest(p, r=r, correction=c("isotropic"))$iso[2]
    for (j in 1:500) {
      thinned_p <- rthin(p, thinning_param)
      k <- Kest(thinned_p, r=r, correction=c("isotropic"))
      thinned_k_hats <- c(thinned_k_hats, k$iso[2])
      normalised_thinned_k_hats <- c(normalised_thinned_k_hats, sqrt(npoints(p))*(k$iso[2]-k_hat))
    }
  }
  return(normalised_thinned_k_hats)
}

r <- c(0.0,0.1)

thinned_k_5_1 <- k_thinned_estimates(0.5)
thinned_k_8_1 <- k_thinned_estimates(0.8)
thinned_k_2_1 <- k_thinned_estimates(0.2)

thinned_k_5_2 <- k_thinned_estimates(0.5)
thinned_k_8_2 <- k_thinned_estimates(0.8)
thinned_k_2_2 <- k_thinned_estimates(0.2)

thinned_k_5_3 <- k_thinned_estimates(0.5)
thinned_k_8_3 <- k_thinned_estimates(0.8)
thinned_k_2_3 <- k_thinned_estimates(0.2)


# par(mfrow=c(1,3))
# qqnorm(thinned_k_2, main="",cex.main=1.5, cex.lab=1.5)
# qqline(thinned_k_2,col=2)
# title("p=0.2",line=1,cex.main=1.5)
# qqnorm(thinned_k_5, main="",cex.main=1.5, cex.lab=1.5)
# qqline(thinned_k_5,col=2)
# title("p=0.5",line=1,cex.main=1.5)
# qqnorm(thinned_k_8, main="",cex.main=1.5, cex.lab=1.5)
# qqline(thinned_k_8,col=2)
# title("p=0.8",line=1,cex.main=1.5)
# title(TeX("\\textbf{QQ-plots of Thinned Estimates, $\\sqrt{n}(\\hat{K}_s(h)$-\\hat{K}(h))}"),outer=TRUE,cex.main=1.7,line=-1)

par(mfrow=c(3,3))
op <- par(mfrow = c(3,3),
          oma = c(5,4,2,1) + 0.1,
          mar = c(2,2,2,2) + 0.1)
qqnorm(thinned_k_2_1, main="",cex.main=1.5, cex.lab=1.5, ylab="Trial 1", xlab="")
qqline(thinned_k_2_1,col=2)
title("Trial 1, p=0.2",line=1,cex.main=1.2)
# title(ylab="Trial 1",line=1.5,outer=TRUE)
qqnorm(thinned_k_5_1, main="",cex.main=1.5, cex.lab=1.5, ylab="", xlab="")
qqline(thinned_k_5_1,col=2)
title("Trial 1, p=0.5",line=1,cex.main=1.2)
qqnorm(thinned_k_8_1, main="",cex.main=1.5, cex.lab=1.5, ylab="", xlab="")
qqline(thinned_k_8_1,col=2)
title("Trial 1, p=0.8",line=1,cex.main=1.2)

qqnorm(thinned_k_2_2, main="",cex.main=1.5, cex.lab=1.5, ylab="Trial 2", xlab="")
qqline(thinned_k_2_2,col=2)
title("Trial 2, p=0.2",line=1,cex.main=1.2)
qqnorm(thinned_k_5_2, main="",cex.main=1.5, cex.lab=1.5, ylab="", xlab="")
qqline(thinned_k_5_2,col=2)
title("Trial 2, p=0.5",line=1,cex.main=1.2)
qqnorm(thinned_k_8_2, main="",cex.main=1.5, cex.lab=1.5, ylab="", xlab="")
qqline(thinned_k_8_2,col=2)
title("Trial 2, p=0.8",line=1,cex.main=1.2)

qqnorm(thinned_k_2_3, main="",cex.main=1.5, cex.lab=1.5, ylab="Trial 3", xlab="")
qqline(thinned_k_2_3,col=2)
title("Trial 3, p=0.2",line=1,cex.main=1.2)
qqnorm(thinned_k_5_3, main="",cex.main=1.5, cex.lab=1.5, ylab="", xlab="")
qqline(thinned_k_5_3,col=2)
title("Trial 3, p=0.5",line=1,cex.main=1.2)
qqnorm(thinned_k_8_3, main="",cex.main=1.5, cex.lab=1.5, ylab="", xlab="")
qqline(thinned_k_8_3,col=2)
title("Trial 3, p=0.8",line=1,cex.main=1.2, ylab="", xlab="")

title(TeX("\\textbf{QQ-plots of Thinned Estimates, $\\sqrt{n}(\\hat{K}_s(h)$-\\hat{K}(h))}"),outer=TRUE,cex.main=1.6,line=1,xlab="Theoretical Quantiles",cex.lab=1.5)
title(ylab="Sample Quantiles",line=1,outer=TRUE,cex.lab=1.5)

