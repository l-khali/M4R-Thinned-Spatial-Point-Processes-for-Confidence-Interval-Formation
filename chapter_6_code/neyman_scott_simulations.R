source("neymann_scott_simulation_function.R")
#source("method_1_version_2.R")
source("method_4_version_2.R")
source("method_3_version_5.R")

ns_tiling_4 <- ns_simulation_tiling(1000, 250, 4, 0.05)
plot(ns_tiling_4[-1,], type="l", lty=3, ylim=c(0.5,1), main="Clustered: Tiling")
ns_tiling_16 <- ns_simulation_tiling(1000, 250, 16, 0.05)
lines(ns_tiling_16[-1,], type="l", lty=2, ylim=c(0.5,1))
ns_tiling_64 <- ns_simulation_tiling(1000, 250, 64, 0.05)
lines(ns_tiling_64[-1,], type="l", lty=1, ylim=c(0.5,1))

ns_approx_4 <- ns_simulation_approximation(1000, 250, 4, 0.05)
plot(ns_approx_4[-1,], type="l", lty=3, ylim=c(0.5,1), main="Clustered: Splitting")
ns_approx_16 <- ns_simulation_approximation(1000, 250, 16, 0.05)
lines(ns_approx_16[-1,], type="l", lty=2, ylim=c(0.5,1))
ns_approx_64 <- ns_simulation_approximation(1000, 250, 64, 0.05)
lines(ns_approx_64[-1,], type="l", lty=1, ylim=c(0.5,1))

ns_subsets_4 <- ns_simulation_subsets5(1000, 250, 4, 0.05)
plot(ns_subsets_4[-1,], type="l", lty=3, ylim=c(0.5,1), main="Clustered: Subsets")
ns_subsets_16 <- ns_simulation_subsets5(1000, 250, 16, 0.05)
lines(ns_subsets_16[-1,], type="l", lty=2, ylim=c(0.5,1))
ns_subsets_64 <- ns_simulation_subsets5(1000, 250, 64, 0.05)
lines(ns_subsets_64[-1,], type="l", lty=1, ylim=c(0.5,1))

ns_marked_4 <- ns_simulation_marked2(50, 250, 4, 0.05)
plot(ns_marked_4[-1,], type="l", lty=3, ylim=c(0.5,1), main="Clustered: Marked point")
ns_marked_16 <- ns_simulation_marked2(50, 250, 16, 0.05)
lines(ns_marked_16[-1,], type="l", lty=2, ylim=c(0.5,1))
ns_marked_64 <- ns_simulation_marked2(50, 250, 64, 0.05)
lines(ns_marked_64[-1,], type="l", lty=1, ylim=c(0.5,1))

ns_thinning_09 <- ns_simulation_thinning(500, 250,0.9,0.05)
ns_thinning_08 <- ns_simulation_thinning(1000, 250,0.9,0.05)
ns_thinning_07 <- ns_simulation_thinning(500, 250,0.7,0.05)
ns_thinning_06 <- ns_simulation_thinning(500, 250,0.6,0.05)
ns_thinning_05 <- ns_simulation_thinning(1000, 250,0.5,0.05)
ns_thinning_04 <- ns_simulation_thinning(500, 250,0.4,0.05)
ns_thinning_03 <- ns_simulation_thinning(500, 250,0.3,0.05)
ns_thinning_02 <- ns_simulation_thinning(1000, 250,0.2,0.05)
ns_thinning_01 <- ns_simulation_thinning(500, 250,0.1,0.05)

par(mgp=c(2,1,0))  
plot(ns_thinning_09[-1,], type="l", lty=1, ylim=c(0.5,1), col=1, xlab="Radius", ylab="Confidence Interval Cover")
lines(ns_thinning_08[-1,], type="l", lty=1, ylim=c(0.5,1), col=2)
lines(ns_thinning_07[-1,], type="l", lty=1, ylim=c(0.5,1), col=3)
lines(ns_thinning_06[-1,], type="l", lty=1, ylim=c(0.5,1), col=4)
lines(ns_thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=5)
lines(ns_thinning_04[-1,], type="l", lty=1, ylim=c(0.5,1), col=6)
lines(ns_thinning_03[-1,], type="l", lty=1, ylim=c(0.5,1), col=7)
lines(ns_thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=8)
lines(ns_thinning_01[-1,], type="l", lty=1, ylim=c(0.5,1), col=9)
legend(0.11, 0.77, legend=c("p=0.9", "p=0.8", "p=0.7", "p=0.6", "p=0.5", "p=0.4", "p=0.3", "p=0.2", "p=0.1"),
       col=c(1,2,3,4,5,6,7,8,9), lty =1, cex=0.8)
abline(h=0.95,lty=2)
title("Cover on Clustered Process", line = 0.3)

 
par(mgp=c(2,1,0))
plot(ns_thinning_08[-1,], type="l", lty=1, ylim=c(0,1), col=1, xlab="Radius", ylab="Confidence Interval Cover")
lines(ns_thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=2)
lines(ns_thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=3)

ns_thinning_08_2 <- ns_thinning_08[-1,]
ns_thinning_05_2 <- ns_thinning_05[-1,]
ns_thinning_02_2 <- ns_thinning_02[-1,]

save(ns_thinning_08_2,file="ns_thinning_08_naive")
save(ns_thinning_05_2,file="ns_thinning_05_naive")
save(ns_thinning_02_2,file="ns_thinning_02_naive")
