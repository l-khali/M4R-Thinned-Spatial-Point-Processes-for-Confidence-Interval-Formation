library(spatstat)
source("soft_core_simulation_function.R")
source("method_2_bootstrapping_spatial_point_patterns.R")

sc_tiling_4_1000 <- sc_simulation_tiling(1000, 250, 4, 0.05)
plot(sc_tiling_4_1000[-1,], type="l", lty=3, ylim=c(0.5,1), main="Soft core: tiling")
sc_tiling_16_1000 <- sc_simulation_tiling(1000, 250, 16, 0.05)
lines(sc_tiling_16_1000[-1,], type="l", lty=2, ylim=c(0.5,1))
sc_tiling_64_1000 <- sc_simulation_tiling(1000, 250, 64, 0.05)
lines(sc_tiling_64_1000[-1,], type="l", lty=1, ylim=c(0.5,1))

sc_approx_4 <- sc_simulation_approximation(1000, 250, 4, 0.05)
plot(sc_approx_4[-1,], type="l", lty=3, ylim=c(0.5,1), main="Unclustered: Splitting")
sc_approx_16 <- sc_simulation_approximation(1000, 250, 16, 0.05)
lines(sc_approx_16[-1,], type="l", lty=2, ylim=c(0.5,1))
sc_approx_64 <- sc_simulation_approximation(1000, 250, 64, 0.05)
lines(sc_approx_64[-1,], type="l", lty=1, ylim=c(0.5,1))

sc_subsets_4_100 <- sc_simulation_subsets(100, 250, 4, 0.05)
plot(sc_subsets_4_100[-1,], type="l", lty=3, ylim=c(0.5,1), main="Soft core: subsets")
sc_tiling_16_100 <- sc_simulation_subsets(100, 250, 16, 0.05)
lines(sc_subsets_16_100[-1,], type="l", lty=2, ylim=c(0.5,1))
sc_subsets_64_100 <- sc_simulation_subsets(100, 250, 64, 0.05)
lines(sc_subsets_64_100[-1,], type="l", lty=1, ylim=c(0.5,1))

sc_marked_4_100 <- sc_simulation_marked(100, 250, 4, 0.05)
plot(sc_marked_4_100[-1,], type="l", lty=3, ylim=c(0.5,1), main="Soft core: marked pt")
sc_marked_16_100 <- sc_simulation_marked(100, 250, 16, 0.05)
lines(sc_marked_16_100[-1,], type="l", lty=2, ylim=c(0.5,1))
sc_marked_64_100 <- sc_simulation_marked(100, 250, 64, 0.05)
lines(sc_marked_64_100[-1,], type="l", lty=1, ylim=c(0.5,1))


sc_thinning_09 <- sc_simulation_thinning(500, 250,0.9,0.05)
sc_thinning_08 <- sc_simulation_thinning(1000, 250,0.9,0.05)
sc_thinning_07 <- sc_simulation_thinning(500, 250,0.7,0.05)
sc_thinning_06 <- sc_simulation_thinning(500, 250,0.6,0.05)
sc_thinning_05 <- sc_simulation_thinning(1000, 250,0.5,0.05)
sc_thinning_04 <- sc_simulation_thinning(500, 250,0.4,0.05)
sc_thinning_03 <- sc_simulation_thinning(500, 250,0.3,0.05)
sc_thinning_02 <- sc_simulation_thinning(1000, 250,0.2,0.05)
sc_thinning_01 <- sc_simulation_thinning(500, 250,0.1,0.05)

par(mgp=c(2,1,0))  
plot(sc_thinning_09[-1,], type="l", lty=1, ylim=c(0.5,1), col=1, xlab="Radius", ylab="Confidence Interval Cover")
lines(sc_thinning_08[-1,], type="l", lty=1, ylim=c(0.5,1), col=2)
lines(sc_thinning_07[-1,], type="l", lty=1, ylim=c(0.5,1), col=3)
lines(sc_thinning_06[-1,], type="l", lty=1, ylim=c(0.5,1), col=4)
lines(sc_thinning_05[-1,], type="l", lty=1, ylim=c(0.5,1), col=5)
lines(sc_thinning_04[-1,], type="l", lty=1, ylim=c(0.5,1), col=6)
lines(sc_thinning_03[-1,], type="l", lty=1, ylim=c(0.5,1), col=7)
lines(sc_thinning_02[-1,], type="l", lty=1, ylim=c(0.5,1), col=8)
lines(sc_thinning_01[-1,], type="l", lty=1, ylim=c(0.5,1), col=9)
legend(0.11, 0.77, legend=c("p=0.9", "p=0.8", "p=0.7", "p=0.6", "p=0.5", "p=0.4", "p=0.3", "p=0.2", "p=0.1"),
       col=c(1,2,3,4,5,6,7,8,9), lty =1, cex=0.8)
abline(h=0.95,lty=2)
title("Cover on Unclustered Process", line = 0.3)

sc_thinning_08_2 <- sc_thinning_08[-1,]
sc_thinning_05_2 <- sc_thinning_05[-1,]
sc_thinning_02_2 <- sc_thinning_02[-1,]

save(sc_thinning_08_2,file="sc_thinning_08_naive")
save(sc_thinning_05_2,file="sc_thinning_05_naive")
save(sc_thinning_02_2,file="sc_thinning_02_naive")


