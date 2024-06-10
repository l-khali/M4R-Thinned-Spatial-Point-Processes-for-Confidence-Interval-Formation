par(mfrow=c(1,3))

par(mgp=c(2,1,0))
plot(thinning_08[-1,], type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, r", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
legend(0.01,0.2,c("p=0.8","p=0.5","p=0.2"),col=c(7,15,22),lwd=c(1.5,1.5,1.5),cex=1.3)
title("Homogenous Poisson",line=1,cex.main=1.5)

par(mgp=c(2,1,0))
plot(ns_thinning_08[-1,], type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, r", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(ns_thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(ns_thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
title("Neyman-Scott",line=1,cex.main=1.5)

par(mgp=c(2,1,0))
plot(sc_thinning_08_2, type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, r", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(sc_thinning_05_2, type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(sc_thinning_02_2, type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
title("Soft Core",line=1,cex.main=1.5)

title("Basic Bootstrap Cover for Ripley's K",line=-1.5,cex.main=1.5,outer=TRUE)



par(mfrow=c(1,3))

par(mgp=c(2,1,0))
plot(thinning_sv_08[-1,], type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, r", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(thinning_sv_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(thinning_sv_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
legend(0.01,0.2,c("p=0.8","p=0.5","p=0.2"),col=c(7,15,22),lwd=c(1.5,1.5,1.5),cex=1.3)
title("Homogenous Poisson",line=1,cex.main=1.5)

par(mgp=c(2,1,0))
plot(ns_thinning_sv_08[-1,], type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, r", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(ns_thinning_sv_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(ns_thinning_sv_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
title("Neyman-Scott",line=1,cex.main=1.5)

par(mgp=c(2,1,0))
plot(sc_thinning_sv_08[-1,], type="l", lty=1, ylim=c(0,1), col=7, lwd=1.5, xlab="Radius, r", ylab="Confidence Interval Cover",cex.lab=1.2)
lines(sc_thinning_sv_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=15, lwd=1.5)
lines(sc_thinning_sv_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=22, lwd=1.5)
abline(h=0.95, col=18, lwd=1.5, lty=2)
title("Soft Core",line=1,cex.main=1.5)

title("Studentized Cover for Ripley's K",line=-1.5,cex.main=1.5,outer=TRUE)

