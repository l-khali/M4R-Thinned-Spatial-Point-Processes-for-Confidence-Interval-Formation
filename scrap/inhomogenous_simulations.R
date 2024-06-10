# wplitting method

inhom_approximation_4 <- poisson_simulation_approximation3(1000,250,4,0.05, inhomogenous=TRUE)
plot(inhom_approximation_4[-1,], ylim = c(0,1), main="Inhomogenous Poisson: splitting", type="l", lty=3)

inhom_approximation_16 <- poisson_simulation_approximation3(1000,250,16,0.05, inhomogenous=TRUE)
lines(inhom_approximation_16[-1,], ylim = c(0,1), type="l", lty=2)

inhom_approximation_64 <- poisson_simulation_approximation3(1000,250,64,0.05, inhomogenous=TRUE)
lines(inhom_approximation_64[-1,], ylim = c(0,1), main="Approximation, N=64")

legend(0.01, 0.72, legend=c("0.1", " 0.2", "0.3","0.4","0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(9,8,7,6,5,4,3,2,1), lty =1, cex=0.7)

# tiling
inhom_tiling_4 <- poisson_simulation_tiling(100, 250, 4, 0.05, inhomogenous = TRUE)
plot(inhom_tiling_4[-1,], ylim = c(0,1), main="Inhomogenous Poisson: tiling", type="l", lty=3)

inhom_tiling_16 <- poisson_simulation_tiling(100, 250, 16, 0.05, inhomogenous = TRUE)
lines(inhom_tiling_16[-1,], ylim = c(0,1), type="l", lty=2)

inhom_tiling_64 <- poisson_simulation_tiling(100, 250, 64, 0.05, inhomogenous = TRUE)
lines(inhom_tiling_64[-1,], ylim = c(0,1), type="l", lty=1)

# subsets
inhom_subsets_4 <- poisson_simulation_subsets5(100,250,4,0.05, inhomogenous=TRUE)
inhom_subsets_16 <- poisson_simulation_subsets5(100,250,16,0.05, inhomogenous=TRUE)
inhom_subsets_64 <- poisson_simulation_subsets5(100,250,64,0.05, inhomogenous=TRUE)
plot(inhom_subsets_4[-1,], ylim=c(0,1), type="l", lty=3, main="Inhomogenous Poisson: subsets")
lines(inhom_subsets_16[-1,], ylim=c(0,1), type="l", lty=2)
lines(inhom_subsets_64[-1,], ylim=c(0,1), type="l", lty=1)

# marked point method
inhom_mrked_4 <- poisson_simulation_marked_point2(100, 250, 4, 0.05, inhomogenous=TRUE)
plot(inhom_mrked_4[-1,], ylim=c(0,1), type="l", lty=3, main="Inhomogenous Poisson: marked point")
inhom_mrked_16 <- poisson_simulation_marked_point2(100, 250, 16, 0.05, inhomogenous=TRUE)
lines(inhom_mrked_16[-1,], ylim=c(0,1), type="l", lty=2)
inhom_mrked_64 <- poisson_simulation_marked_point2(100, 250, 64, 0.05, inhomogenous=TRUE)
lines(inhom_mrked_64[-1,], ylim=c(0,1), type="l", lty=1)