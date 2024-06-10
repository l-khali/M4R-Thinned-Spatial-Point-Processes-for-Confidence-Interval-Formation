library(spatstat)
library(hash)

K_var_inhom <- function(intensity, N = 4, B=1000, radius = seq(0.0,0.14,0.01), thinning_param = 1, window = owin(c(0,1),c(0,1))) {
  # vars <- data.frame(radius)
  p <- rpoispp(intensity, win=window)
  if (thinning_param < 1) {
    unif <- runif(npoints(p), 0, 1)
    p_df <- as.data.frame(p)
    subprocess_df <- p_df[which(unif < thinning_param),]
    rownames(subprocess_df) <- seq(1, nrow(subprocess_df))
    p <- as.ppp(subprocess_df, window)
  }
  npoints <- npoints(p)
  r <- t(as.data.frame(replicate(npoints,radius)))
  weights <- edge.Ripley(p,pairdist(p))
  a <- 1/N
  increment <- 1/sqrt(N)
  window <- owin(c(0, increment), c(0,increment))
  k_vals <- data.frame(radius)
    
  marks <- as.data.frame(matrix(0, ncol = length(radius), nrow = npoints))
  datadf <- as.data.frame(p)
    
  for (point in 1:(nrow(datadf)-1)) {
    x1 <- datadf[point,1]
    y1 <- datadf[point,2]
    for (other_point in (point+1):(nrow(datadf))) {
      x2 <- datadf[other_point,1]
      y2 <- datadf[other_point,2]
      for (radius_index in 1:15) {
        current_radius <- seq(0.0, 0.14, 0.01)[radius_index]
        if ((x1-x2)**2+(y1-y2)**2 <= current_radius**2) {
          marks[point,radius_index:15] <- marks[point,radius_index:15] + rep(weights[point, other_point], each=length(radius_index:15))
          marks[other_point,radius_index:15] <- marks[other_point,radius_index:15] + rep(weights[other_point, point], each=length(radius_index:15))
          break
        }
      }
    }
  }
  marks(p) <- marks
    
  # only for N=4 for now
  blocks = hash()
  block_count <- 0
  for (xstart in head(seq(0,1,length.out = sqrt(N) + 1),-1)) {
    for (ystart in head(seq(0,1,length.out = sqrt(N) + 1),-1)) {
      block_count = block_count +  1
      blocks[[as.character(block_count)]] = subset(p, xstart <= x & x < xstart + increment & ystart <= y & y < ystart + increment)
    }
  }
  
  bootstrapped_K <- data.frame(radius)
  for (bootstrap in 1:B) {
    current_bootstrap_K <- c()
    sampled_blocks <- sample(c(1,2,3,4), 4, replace=TRUE)
    current_K <- data.frame(radius)
    
    for (block in sampled_blocks) {
      block_k <- c()
      for (r in 1:length(radius)) {
        block_k <- c(block_k,sum(marks(blocks[[as.character(block)]])[,r]))
      }
      current_K <- cbind(current_K, block_k)
    }
    
    for (r in 1:length(radius)) {
      current_bootstrap_K <- c(current_bootstrap_K, sum(current_K[r,-1]))
    }
    bootstrapped_K <- cbind(bootstrapped_K, current_bootstrap_K)
  }
  variances <- c()
  for (r in 1:length(radius) ) {
    variances <- c(variances, var(as.numeric(bootstrapped_K[r,-1])))
  }
  
  return(variances) 
}
