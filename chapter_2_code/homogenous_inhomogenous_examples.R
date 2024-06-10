library(spatstat)

hom <- rpoispp(350)
inhom <- rpoispp(function(x,y){500*(x^2+y^2)})
inhom2 <- rpoispp(function(x,y){15*exp(5*x)})

plot(hom)
plot(inhom)
plot(inhom2)

par(mfrow=c(1,3), mar = c(1, 1, 1, 1))
plot(hom,main="")
title("Homogenous",cex.main=1.5,line=-1)
plot(inhom,main="")
title("Inhomogenous",cex.main=1.5,line=-1)
plot(inhom2,main="")
title("Inhomogenous",cex.main=1.5,line=-1)

