par(mfrow=c(1,3))

plot(thinning_08_2,ylim=c(0.6,1),type="l",col=7,lwd=1.5,xlab=TeX('Radius, $h$'),ylab="Confidence Interval Cover")
lines(thinning_05_2,ylim=c(0.5,1),col=15,lwd=1.5)
lines(thinning_02_2,ylim=c(0.5,1),col=22,lwd=1.5)
abline(h=0.95, col=18,lwd=1.5,lty=2)
legend(0.02,0.7,c("p=0.8","p=0.5","p=0.2"),col=c(7,15,22),lty=c(1,1,1),lwd=c(1.5,1.5,1.5),cex=0.9)
title("Homogenous Poisson",line=0.4)

plot(ns_thinning_08_2,ylim=c(0.2,1),type="l",col=7,lwd=1.5,xlab=TeX('Radius, $h$'),ylab="Confidence Interval Cover")
lines(ns_thinning_05_2,ylim=c(0.5,1),col=15,lwd=1.5)
lines(ns_thinning_02_2,ylim=c(0.5,1),col=22,lwd=1.5)
abline(h=0.95, col=18,lwd=1.5,lty=2)
title("Neyman-Scott",line=0.4)

plot(sc_thinning_08_2,ylim=c(0,1),type="l",col=7,lwd=1.5,xlab=TeX('Radius, $h$'),ylab="Confidence Interval Cover")
lines(sc_thinning_05_2,ylim=c(0.5,1),col=15,lwd=1.5)
lines(sc_thinning_02_2,ylim=c(0.5,1),col=22,lwd=1.5)
abline(h=0.95, col=18,lwd=1.5,lty=2)
title("Soft Core",line=0.4)

title("Basic Boostrap CI Cover: Spatial Processes", line = -2, outer = TRUE,cex.main=1.5)

