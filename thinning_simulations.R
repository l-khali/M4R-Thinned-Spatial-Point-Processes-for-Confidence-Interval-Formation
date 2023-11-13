# homogenous poisson simulations

thinning_sv_09 <- poisson_simulation_thinning_sv(100, 250, 0.9, 0.05)
thinning_sv_08 <- poisson_simulation_thinning_sv(100, 250, 0.8, 0.05)
thinning_sv_07 <- poisson_simulation_thinning_sv(100, 250, 0.7, 0.05)
thinning_sv_06 <- poisson_simulation_thinning_sv(100, 250, 0.6, 0.05)
thinning_sv_05 <- poisson_simulation_thinning_sv(100, 250, 0.5, 0.05)
thinning_sv_04 <- poisson_simulation_thinning_sv(100, 250, 0.4, 0.05)
thinning_sv_03 <- poisson_simulation_thinning_sv(100, 250, 0.3, 0.05)
thinning_sv_02 <- poisson_simulation_thinning_sv(100, 250, 0.2, 0.05)
thinning_sv_01 <- poisson_simulation_thinning_sv(100, 250, 0.1, 0.05)

plot(thinning_sv_09[-1,], ylim=c(0.5,1), type="l", col=1, main="Poisson: Thinning")
lines(thinning_sv_08[-1,], type="l", col=2)
lines(thinning_sv_07[-1,], type="l", col=3)
lines(thinning_sv_06[-1,], type="l", col=4)
lines(thinning_sv_05[-1,], type="l", col=5)
lines(thinning_sv_04[-1,], type="l", col=6)
lines(thinning_sv_03[-1,], type="l", col=7)
lines(thinning_sv_02[-1,], type="l", col=8)
lines(thinning_sv_01[-1,], type="l", col=9)
legend(0.01, 0.72, legend=c("0.1", " 0.2", "0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(9,8,7,6,5,4,3,2,1), lty =1, cex=0.7)
abline(h = 0.95, col=1, lty=2)

# inhomogenous poisson simulations

thinning_sv_09_inhom <- poisson_simulation_thinning_sv(100, 250, 0.9, 0.05, inhomogenous = TRUE)
thinning_sv_08_inhom <- poisson_simulation_thinning_sv(100, 250, 0.8, 0.05, inhomogenous = TRUE)
thinning_sv_07_inhom <- poisson_simulation_thinning_sv(100, 250, 0.7, 0.05, inhomogenous = TRUE)
thinning_sv_06_inhom <- poisson_simulation_thinning_sv(100, 250, 0.6, 0.05, inhomogenous = TRUE)
thinning_sv_05_inhom <- poisson_simulation_thinning_sv(100, 250, 0.5, 0.05, inhomogenous = TRUE)
thinning_sv_04_inhom <- poisson_simulation_thinning_sv(100, 250, 0.4, 0.05, inhomogenous = TRUE)
thinning_sv_03_inhom <- poisson_simulation_thinning_sv(100, 250, 0.3, 0.05, inhomogenous = TRUE)
thinning_sv_02_inhom <- poisson_simulation_thinning_sv(100, 250, 0.2, 0.05, inhomogenous = TRUE)
thinning_sv_01_inhom <- poisson_simulation_thinning_sv(100, 250, 0.1, 0.05, inhomogenous = TRUE)

plot(thinning_sv_09_inhom[-1,], ylim=c(0,1), type="l", col=1, main="Inhomogenous Poisson: Thinning (using Loh bootstrapping, N=16)")
lines(thinning_sv_08_inhom[-1,], type="l", col=2)
lines(thinning_sv_07_inhom[-1,], type="l", col=3)
lines(thinning_sv_06_inhom[-1,], type="l", col=4)
lines(thinning_sv_05_inhom[-1,], type="l", col=5)
lines(thinning_sv_04_inhom[-1,], type="l", col=6)
lines(thinning_sv_03_inhom[-1,], type="l", col=7)
lines(thinning_sv_02_inhom[-1,], type="l", col=8)
lines(thinning_sv_01_inhom[-1,], type="l", col=9)
legend(0.01, 0.4, legend=c("0.1", " 0.2", "0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(9,8,7,6,5,4,3,2,1), lty =1, cex=0.7)
abline(h = 0.95, col=1, lty=2)

# soft core simulations

sc_thinning_sv_09 <- sc_simulation_thinning_sv(100, 250, 0.9, 0.05)
sc_thinning_sv_08 <- sc_simulation_thinning_sv(100, 250, 0.8, 0.05)
sc_thinning_sv_07 <- sc_simulation_thinning_sv(100, 250, 0.7, 0.05)
sc_thinning_sv_06 <- sc_simulation_thinning_sv(100, 250, 0.6, 0.05)
sc_thinning_sv_05 <- sc_simulation_thinning_sv(100, 250, 0.5, 0.05)
sc_thinning_sv_04 <- sc_simulation_thinning_sv(100, 250, 0.4, 0.05)
sc_thinning_sv_03 <- sc_simulation_thinning_sv(100, 250, 0.3, 0.05)
sc_thinning_sv_02 <- sc_simulation_thinning_sv(100, 250, 0.2, 0.05)
sc_thinning_sv_01 <- sc_simulation_thinning_sv(100, 250, 0.1, 0.05)

plot(sc_thinning_sv_09[-1,], ylim=c(0.5,1), type="l", col=1, main="Unclustered: Thinnning")
lines(sc_thinning_sv_08[-1,], type="l", col=2)
lines(sc_thinning_sv_07[-1,], type="l", col=3)
lines(sc_thinning_sv_06[-1,], type="l", col=4)
lines(sc_thinning_sv_05[-1,], type="l", col=5)
lines(sc_thinning_sv_04[-1,], type="l", col=6)
lines(sc_thinning_sv_03[-1,], type="l", col=7)
lines(sc_thinning_sv_02[-1,], type="l", col=8)
lines(sc_thinning_sv_01[-1,], type="l", col=9)
legend(0.01, 0.72, legend=c("0.1", " 0.2", "0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(9,8,7,6,5,4,3,2,1), lty =1, cex=0.7)
abline(h = 0.95, col=1, lty=2)

# neyman scott simulations

ns_thinning_sv_09 <- ns_simulation_thinning_sv(100, 250, 0.9, 0.05)
ns_thinning_sv_08 <- ns_simulation_thinning_sv(100, 250, 0.8, 0.05)
ns_thinning_sv_07 <- ns_simulation_thinning_sv(100, 250, 0.7, 0.05)
ns_thinning_sv_06 <- ns_simulation_thinning_sv(100, 250, 0.6, 0.05)
ns_thinning_sv_05 <- ns_simulation_thinning_sv(100, 250, 0.5, 0.05)
ns_thinning_sv_04 <- ns_simulation_thinning_sv(100, 250, 0.4, 0.05)
ns_thinning_sv_03 <- ns_simulation_thinning_sv(100, 250, 0.3, 0.05)
ns_thinning_sv_02 <- ns_simulation_thinning_sv(100, 250, 0.2, 0.05)
ns_thinning_sv_01 <- ns_simulation_thinning_sv(100, 250, 0.1, 0.05)

plot(ns_thinning_sv_09[-1,], ylim=c(0,1), type="l", col=1, main="Clustered: Thinning")
lines(ns_thinning_sv_08[-1,], type="l", col=2)
lines(ns_thinning_sv_07[-1,], type="l", col=3)
lines(ns_thinning_sv_06[-1,], type="l", col=4)
lines(ns_thinning_sv_05[-1,], type="l", col=5)
lines(ns_thinning_sv_04[-1,], type="l", col=6)
lines(ns_thinning_sv_03[-1,], type="l", col=7)
lines(ns_thinning_sv_02[-1,], type="l", col=8)
lines(ns_thinning_sv_01[-1,], type="l", col=9)
abline(h = 0.95, col=1, lty=2)
legend(0.11, 1, legend=c("0.1", " 0.2", "0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(9,8,7,6,5,4,3,2,1), lty =1, cex=0.7)


