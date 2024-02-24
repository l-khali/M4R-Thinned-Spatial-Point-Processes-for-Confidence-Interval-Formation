library(spatstat)
library(hash)

K_var <- function(radius, npoints, window = owin(c(0,1),c(0,1))) {
  # This is the variance of K on homogenous data as stated in Lagache, 2013, pg2
  omega <- area(window)
  beta <- pi * (radius^2) / omega
  gamma <- perimeter(window) * radius / omega
  return(((2*(omega^2)*beta) * (1+0.305*gamma+beta*(-1+0.0132*npoints*gamma)))/(npoints^2))
}

variances <- c()
for (npoints in 10:280) {
  variances <- c(variances, K_var(0.01,npoints))
}

scaling_constant <- function(thinning_param, npoints, radius=seq(0,0.14,0.01), mean=FALSE) {
  if (mean==TRUE) {
    return (mean((K_var(radius,round(npoints*thinning_param))/K_var(radius,npoints))[-1]))
  } else {
    return(K_var(radius,round(npoints*thinning_param))/K_var(radius,npoints))
  }
}

# K_var_inhom <- function(intensity, radius = seq(0.0,0.14,0.01), thinning_param = 1, window = owin(c(0,1),c(0,1))) {
#   k_list <- data.frame(r)
#   for (sim in 1:1000) {
#     p <- rpoispp(intensity)
#     if (thinning_param < 1) {
#       unif <- runif(npoints(p), 0, 1)
#       p_df <- as.data.frame(p)
#       subprocess_df <- p_df[which(unif < thinning_param),]
#       p <- as.ppp(subprocess_df, window)
#     }
#     k_list <- cbind(k_list, as.data.frame(Kinhom(p, lambda=intensity, r = radius, correction=c("isotropic")))["iso"])
#   }
#   vars <- c()
#   for (r in 1:length(radius)) {
#     vars <- c(vars, var(as.numeric(k_list[r,-1])))
#   }
#   return(vars)
# }
# 
# K_var_inhom <- function(intensity, N = 4, radius = seq(0.0,0.14,0.01), thinning_param = 1, window = owin(c(0,1),c(0,1))) {
#   vars <- data.frame(radius)
#   for (sim in 1:10) {
#     p <- rpoispp(intensity, win=window)
#     if (thinning_param < 1) {
#       unif <- runif(npoints(p), 0, 1)
#       p_df <- as.data.frame(p)
#       subprocess_df <- p_df[which(unif < thinning_param),]
#       rownames(subprocess_df) <- seq(1, nrow(subprocess_df))
#       p <- as.ppp(subprocess_df, window)
#     }
#     plot(p)
#     npoints <- npoints(p)
#     r <- t(as.data.frame(replicate(npoints,radius)))
#     weights <- edge.Ripley(p,pairdist(p))
#     a <- 1/N
#     increment <- 1/sqrt(N)
#     window <- owin(c(0, increment), c(0,increment))
#     k_vals <- data.frame(radius)
#     
#     marks <- as.data.frame(matrix(0, ncol = length(radius), nrow = npoints))
#     datadf <- as.data.frame(p)
#     
#     for (point in 1:(nrow(datadf)-1)) {
#       x1 <- datadf[point,1]
#       y1 <- datadf[point,2]
#       for (other_point in (point+1):(nrow(datadf))) {
#         x2 <- datadf[other_point,1]
#         y2 <- datadf[other_point,2]
#         for (radius_index in 1:15) {
#           current_radius <- seq(0.0, 0.14, 0.01)[radius_index]
#           if ((x1-x2)**2+(y1-y2)**2 <= current_radius**2) {
#             marks[point,radius_index:15] <- marks[point,radius_index:15] + rep(weights[point, other_point], each=length(radius_index:15))
#             marks[other_point,radius_index:15] <- marks[other_point,radius_index:15] + rep(weights[other_point, point], each=length(radius_index:15))
#             break
#           }
#         }
#       }
#     }
#     marks(p) <- marks
#     
#     # only for N=4 for now
#     blocks = hash()
#     block_count <- 0
#     for (xstart in head(seq(0,1,length.out = sqrt(N) + 1),-1)) {
#       for (ystart in head(seq(0,1,length.out = sqrt(N) + 1),-1)) {
#         block_count = block_count +  1
#         blocks[[as.character(block_count)]] = subset(p, xstart <= x & x < xstart + increment & ystart <= y & y < ystart + increment)
#       }
#     }
#     print(blocks)
#     
#     s1 <- c()
#     s2 <- rep(0,length(radius))
#     s3 <- rep(0,length(radius))
#     s4 <- rep(0,length(radius))
#     s5 <- rep(0,length(radius))
#     
#     for (r in 1:length(radius)) {
#       s1 <- c(s1, sum(marks(p)[,r]^2))
#     }
#     
#     for (block1 in 1:N) {
#       block_marks <- marks(blocks[[as.character(block1)]])
#       rownames(block_marks) <- seq(1,nrow(block_marks))
#         
#       for (r in 1:length(radius)) {
#         s4 <- c(s4, sum(block_marks[,r])^2)
#       }
#         
#       for (point1 in 1:(nrow(block_marks)-1)) {
#         for (point2 in (point1+1):nrow(block_marks)) {
#           s2 <- s2 + block_marks[point1,] * block_marks[point2,]
#         }
#       }
#     }
#     
#     for (block1 in 1:(N-1)) {
#         for (block2 in (block1+1):N) {
#           block_marks1 <- marks(blocks[[as.character(block1)]])
#           rownames(block_marks1) <- seq(1,nrow(block_marks1))
#           block_marks2 <- marks(blocks[[as.character(block2)]])
#           rownames(block_marks2) <- seq(1,nrow(block_marks2))
#           s5 <- s5 + colSums(block_marks1) * colSums(block_marks2)
#           for (point1 in 1:nrow(block_marks1)) {
#             for (point2 in 1:nrow(block_marks2)) {
#               s3 <- s3 + block_marks1[point1,] * block_marks2[point2,]
#             }
#           }
#   
#         }
#         
#     }
#     vars <- cbind(vars, t((s1+s2+s3-s4-s5)/(area(window)^2)))
#   }
#   average_vars <- c()
#   for (r in 1:length(radius)) {
#     average_vars <- c(average_vars, mean(as.numeric(vars[r,-1])))
#   }
#   return(average_vars)
# }

scaling_constant_inhom <- function(thinning_param, intensity, radius=seq(0,0.14,0.01), mean=FALSE) {
  if (mean==TRUE) {
    return(mean(as.numeric((K_var_inhom(intensity, radius=radius, thinning_param = thinning_param)/K_var_inhom(intensity, radius=radius, thinning_param = 1))[-1])))
  } else {
    return((K_var_inhom(intensity, radius=radius, thinning_param = thinning_param)/K_var_inhom(intensity, radius=radius, thinning_param = 1)))
  }
}

# homogenous

plot(seq(0,0.14,0.01), scaling_constant(0.9, 250), type="l", lty=1, col=1, ylim=c(0,100), xlab="radius", ylab="scaling factor", main="radius vs. scaling factor")
lines(seq(0,0.14,0.01), scaling_constant(0.8, 250), lty=1, col=2)
lines(seq(0,0.14,0.01), scaling_constant(0.7, 250), lty=1, col=3)
lines(seq(0,0.14,0.01), scaling_constant(0.6, 250), lty=1, col=4)
lines(seq(0,0.14,0.01), scaling_constant(0.5, 250), lty=1, col=5)
lines(seq(0,0.14,0.01), scaling_constant(0.4, 250), lty=1, col=6)
lines(seq(0,0.14,0.01), scaling_constant(0.3, 250), lty=1, col=7)
lines(seq(0,0.14,0.01), scaling_constant(0.2, 250), lty=1, col=8)
lines(seq(0,0.14,0.01), scaling_constant(0.1, 250), lty=1, col=9)
legend(0, 95, legend=c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(9,8,7,6,5,4,3,2,1), lty =1, cex=0.7)


plot(seq(0,0.14,0.01), scaling_constant(0.9, 250), type="l", lty=1, col=1, ylim=c(0,5), xlab="radius", ylab="scaling factor", main="radius vs. scaling factor (zoomed in)")
lines(seq(0,0.14,0.01), scaling_constant(0.8, 250), lty=1, col=2)
lines(seq(0,0.14,0.01), scaling_constant(0.7, 250), lty=1, col=3)
lines(seq(0,0.14,0.01), scaling_constant(0.6, 250), lty=1, col=4)
lines(seq(0,0.14,0.01), scaling_constant(0.5, 250), lty=1, col=5)
legend(0, 5, legend=c("0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(5,4,3,2,1), lty =1, cex=0.7)

# inhomogenous
plot(seq(0,0.14,0.01), scaling_constant_inhom(0.9, intensity2), type="l", lty=1, col=1, ylim=c(0,10), xlab="radius", ylab="scaling factor", main="radius vs. scaling factor (inhomogenous)")
lines(seq(0,0.14,0.01), scaling_constant_inhom(0.8, intensity2), lty=1, col=2)
lines(seq(0,0.14,0.01), scaling_constant_inhom(0.7, intensity2), lty=1, col=3)
lines(seq(0,0.14,0.01), scaling_constant_inhom(0.6, intensity2), lty=1, col=4)
lines(seq(0,0.14,0.01), scaling_constant_inhom(0.5, intensity2), lty=1, col=5)
lines(seq(0,0.14,0.01), scaling_constant_inhom(0.4, intensity2), lty=1, col=6)
lines(seq(0,0.14,0.01), scaling_constant_inhom(0.3, intensity2), lty=1, col=7)
lines(seq(0,0.14,0.01), scaling_constant_inhom(0.2, intensity2), lty=1, col=8)
lines(seq(0,0.14,0.01), scaling_constant_inhom(0.1, intensity2), lty=1, col=9)
legend(0, 10, legend=c("0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"),
       col=c(9,8,7,6,5,4,3,2,1), lty =1, cex=0.7)

