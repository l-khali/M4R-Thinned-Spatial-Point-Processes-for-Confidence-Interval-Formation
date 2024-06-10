# Creates a plot to demonstrate how the timing method works with toroidal wrapping

# Generate a homogeneous poisson process with lambda=50
data <- rpoispp(50)
plot(data,pch=20, cex=1.5,main="")
title("Original Poisson Process", line=-0.6, cex.main=2)

# Specifying increments used for tiling method when N = 4
new_process <- data.frame()
increments <- data.frame()
for (i in 1:sqrt(4)) {
  for (j in 1:sqrt(4)) {
    increments <- rbind(increments, c((i-1)*0.5, (j-1)*0.5))
  }
}

# Drawing "randomly" placed tile
xstart <- 0.1
ystart <- 0.15
xend <- 0.6
yend <- 0.65

segments(xstart, ystart, xend, ystart, lty=2)
segments(xstart, ystart, xstart, yend, lty=2)
segments(xstart, yend, xend, yend, lty=2)
segments(xend, ystart, xend, yend, lty=2)
text(0.12, 0.63, "a", cex=1.5)

# Data in this tile to predetermined position
if (xend <= 1 & yend <= 1) {
  subregion <- as.data.frame(subset(data, xstart <= x & x < xend & ystart <= y & y < yend))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[1,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[1,2])
  new_process <- rbind(new_process, subregion)
} else if (xend > 1 & yend <= 1) {
  subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & ystart <= y & y < yend))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[1,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[1,2])
  new_process <- rbind(new_process, subregion)
} else if (xend <= 1 & yend > 1) {
  subregion <- as.data.frame(subset(data, xstart <= x & x < xend & (ystart <= y | y < yend %% 1)))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[1,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[1,2])
  new_process <- rbind(new_process, subregion)
} else {
  subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y | y < yend %% 1)))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[1,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[1,2])
  new_process <- rbind(new_process, subregion)
}

# Tiling process repeated three more times

xstart <- 0.75
ystart <- 0.3
xend <- 1.25
yend <- 0.8

segments(0.75,0.3,1,0.3, lty=2, col=2)
segments(0.75,0.8,1,0.8, lty=2, col=2)
segments(0,0.3,0.25,0.3, lty=2, col=2)
segments(0,0.8,0.25,0.8, lty=2, col=2)
segments(0.25,0.3,0.25,0.8, lty=2, col=2)
segments(0.75,0.3,0.75,0.8, lty=2, col=2)
text(0.77, 0.78, "b", cex=1.5)

if (xend <= 1 & yend <= 1) {
  subregion <- as.data.frame(subset(data, xstart <= x & x < xend & ystart <= y & y < yend))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[2,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[2,2])
  new_process <- rbind(new_process, subregion)
} else if (xend > 1 & yend <= 1) {
  subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & ystart <= y & y < yend))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[2,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[2,2])
  new_process <- rbind(new_process, subregion)
} else if (xend <= 1 & yend > 1) {
  subregion <- as.data.frame(subset(data, xstart <= x & x < xend & (ystart <= y | y < yend %% 1)))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[2,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[2,2])
  new_process <- rbind(new_process, subregion)
} else {
  subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y | y < yend %% 1)))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[2,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[2,2])
  new_process <- rbind(new_process, subregion)
}


xstart <- 0.35
ystart <- 0.7
xend <- 0.85
yend <- 1.2

segments(0.35, 0.7, 0.85, 0.7, lty=2, col=2)
segments(0.35, 0.2, 0.85, 0.2, lty=2, col=2)
segments(0.35,0.7,0.35,1, lty=2, col=2)
segments(0.85,0.7,0.85,1, lty=2, col=2)
segments(0.35,0,0.35,0.2, lty=2, col=2)
segments(0.85,0,0.85,0.2, lty=2, col=2)
text(0.82, 0.18, "c", cex=1.5)

if (xend <= 1 & yend <= 1) {
  subregion <- as.data.frame(subset(data, xstart <= x & x < xend & ystart <= y & y < yend))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[3,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[3,2])
  new_process <- rbind(new_process, subregion)
} else if (xend > 1 & yend <= 1) {
  subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & ystart <= y & y < yend))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[3,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[3,2])
  new_process <- rbind(new_process, subregion)
} else if (xend <= 1 & yend > 1) {
  subregion <- as.data.frame(subset(data, xstart <= x & x < xend & (ystart <= y | y < yend %% 1)))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[3,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[3,2])
  new_process <- rbind(new_process, subregion)
} else {
  subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y | y < yend %% 1)))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[3,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[3,2])
  new_process <- rbind(new_process, subregion)
}

xstart <- 0.48
ystart <- 0.1
xend <- 0.98
yend <- 0.6

segments(xstart, ystart, xend, ystart, lty=2)
segments(xstart, ystart, xstart, yend, lty=2)
segments(xstart, yend, xend, yend, lty=2)
segments(xend, ystart, xend, yend, lty=2)
text(0.5, 0.58, "d", cex=1.5)

if (xend <= 1 & yend <= 1) {
  subregion <- as.data.frame(subset(data, xstart <= x & x < xend & ystart <= y & y < yend))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[4,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[4,2])
  new_process <- rbind(new_process, subregion)
} else if (xend > 1 & yend <= 1) {
  subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & ystart <= y & y < yend))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[4,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[4,2])
  new_process <- rbind(new_process, subregion)
} else if (xend <= 1 & yend > 1) {
  subregion <- as.data.frame(subset(data, xstart <= x & x < xend & (ystart <= y | y < yend %% 1)))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[4,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[4,2])
  new_process <- rbind(new_process, subregion)
} else {
  subregion <- as.data.frame(subset(data, (xstart <= x | x < xend %% 1) & (ystart <= y | y < yend %% 1)))
  subregion["x"] <- lapply(subregion["x"], function(x) ((x - xstart) %% 1) + increments[4,1])
  subregion["y"] <- lapply(subregion["y"], function(y) ((y - ystart) %% 1) + increments[4,2])
  new_process <- rbind(new_process, subregion)
}

window <- owin(c(0, 1), c(0,1))
new_ppp <- as.ppp(new_process, window)

plot(new_ppp,pch=20, cex=1.5, main="")
title("Bootstrapped Process (Tiling)", line=-0.6, cex.main=2)
segments(0.5,0,0.5,1,lty=2)
segments(0,0.5,1,0.5,lty=2)
text(0.02, 0.48, "a", cex=1.5)
text(0.02, 0.98, "b", cex=1.5)
text(0.52, 0.48, "c", cex=1.5)
text(0.52, 0.98, "d", cex=1.5)

