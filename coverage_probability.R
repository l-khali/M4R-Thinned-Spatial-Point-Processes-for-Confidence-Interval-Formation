library(spatstat)

coverage_probability <- function(confidences, actual) {
  covered <- 0
  print(nrow(confidences))
  for (i in 1:nrow(confidences)) {
    if (confidences[i,1] <= actual & actual <= confidences[i,2]) {
      covered <- covered + 1
    }
  }
  return(covered/nrow(confidences))
}
