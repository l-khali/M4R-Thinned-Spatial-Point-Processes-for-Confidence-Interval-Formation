library(spatstat)
library(abind)
library(pals)
library(latex2exp)

intensity_func <- function(x,y) {
  return(exp(a+b*x))
}

kernel_cover <- function(nsim, a, b, alpha, thinning_param, B) {
  
  cover <- matrix(rep(0,128*128),128,128)
  scalar <- thinning_param/(1-thinning_param)
  t <- qt((1-alpha/2), df = B-1)
  
  ci_func <- function(v) {
    return(t*sqrt(v*scalar))
  }
  
  for (sim in 1:nsim) {
    print(paste("Current simulation:",sim))
    p <- rpoispp(intensity_func,win=square(1))
    d_hat <- density(p,diggle=TRUE,at="pixels",positive = TRUE)
    xs <- d_hat$xcol
    ys <- d_hat$yrow
    for (bsim in 1:B) {
      subp <- rthin(p, thinning_param)
      d <- density(subp,diggle=TRUE,at="pixels",positive = TRUE)/thinning_param
      if (bsim == 1) {
        ds <- as.matrix(d)
      } else {
        ds <- abind(ds, as.matrix(d), along=3)
      }
    }
    kernel_vars <- apply(ds, c(1,2), var)
    ci <- ci_func(kernel_vars)
    upper <- d_hat - mean(d_hat) + ci
    lower <- d_hat - mean(d_hat) - ci
    
    true_kernel <- t(outer(ys,xs,intensity_func))
    true_kernel <- true_kernel - mean(true_kernel)
    cover <- cover + (as.matrix(true_kernel < upper & true_kernel > lower))/nsim
    
  }
  return(cover)
}


kernel_cover_8 <- kernel_cover(200,3,3,0.05,0.8,500)
kernel_cover_5 <- kernel_cover(200,3,3,0.05,0.5,500)
kernel_cover_2 <- kernel_cover(200,3,3,0.05,0.2,500)




cm1 <- colourmap(cubicyf(100), range=c(0,1))
cm2 <- colourmap(cubicyf(100), range=c(0.93,0.97))

par(mfrow = c(3, 2))
op <- par(mfrow = c(3,2),
          oma = c(0.5,2,4,0.5) + 0.1,
          mar = c(0,0,0,2) + 0.1,
          mai = c(0.3,0,0,0.5))

plot(im(kernel_cover_2),col=cm1,ncolours=30, main="")
txt <- paste0("$p = 0.2$")
title(TeX(txt), line=0.5, cex.main=1.5,outer=TRUE)
plot(im(kernel_cover_2),col=cm2,ncolours=30, main="")

plot(im(kernel_cover_5),col=cm1,ncolours=30, main="")
txt <- paste0("$p = 0.5$")
title(TeX(txt), line=-14.5, cex.main=1.5,outer=TRUE)
plot(im(kernel_cover_5),col=cm2,ncolours=30, main="")

plot(im(kernel_cover_8),col=cm1,ncolours=30, main="")
txt <- paste0("$p = 0.8$")
title(TeX(txt), line=-29.5, cex.main=1.5,outer=TRUE)
plot(im(kernel_cover_8),col=cm2,ncolours=30, main="")

title(TeX("Cover of Kernel Estimate Using Thinning"),line=2.5,outer=TRUE, cex.main=1.6)
title(xlab="x", outer=TRUE, line=-1.5, cex.lab=1.5)
title(ylab="y", outer=TRUE, line=1, cex.lab=1.5)


par(mfrow=c(1,3), mar = c(1, 1, 3, 1))
cm3 <- colourmap(cubicyf(100), range=c(0,400))
a <- 3
b <- 3
p <- rpoispp(intensity_func,win=square(1))
d_hat <- density(p,diggle=TRUE,at="pixels",positive = TRUE)
xs <- d_hat$xcol
ys <- d_hat$yrow
true_kernel <- t(outer(ys,xs,intensity_func))
plot(im(true_kernel), col=cm3, main="")
title("True Intensity",line=0,cex.main=1.4)
plot(d_hat, col=cm3, main="")
title("Kernel Estimate",line=0,cex.main=1.4)
plot(abs(true_kernel-d_hat)/true_kernel, col=brewer.ylorrd(100), main="")
title("Relative \n Absolute Error",line=0,cex.main=1.4)





