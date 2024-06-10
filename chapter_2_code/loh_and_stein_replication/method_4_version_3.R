library(spatstat)
library(ecespa)

marked_point_method3 <- function(data, N, alpha, R=99) {
  r <- seq(0.0, 0.14, 0.01)
  weights <- edge.Ripley(data, pairdist(data))
  marked_process <- data
  marks <- as.data.frame(matrix(0, ncol = length(r), nrow = npoints(data)))
  # populating the marks df is very very very slow
  for (point in 1:npoints(data)) {
    for (radius in 1:length(r)) {
      data$marks <- weights[point,]
      marks[point, radius] <- marksum(data, r[radius])$marksum[point]
    }
  }
  return(marks)
}