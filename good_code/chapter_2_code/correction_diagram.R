library(spatstat)
library(plotrix)
library(graphics)
library(pals)

pois <- rpoispp(250)
plot(pois)

matern <- rMatClust(10,0.1,25)
plot(matern)

thomas <- rThomas(10,0.1,25)
plot(thomas)

sc <- simulate_soft_core()
plot(sc)

par(mfrow=c(1,4),mar=c(1,0.25,0,0.25))
plot(pois,main="")
title("Homogenous Poisson", line=-18)
plot(matern,main="")
title("Matérn Cluster Process", line=-18)
plot(thomas,main="")
title("Thomas Process", line=-18)
plot(sc,main="")
title("Soft Core Process", line=-18)


pois2 <- rpoispp(100,win=square(2))

palette(alphabet(26))
plot(subset(pois2,owin(c(0.25,1.75),c(0.25,1.75))),show.window	=FALSE,main="")
rect(0.5, 0.5, 1.5, 1.5, density = NULL, angle = 45,
     col = NA, border = NULL, lty = par("lty"), lwd = 2)
text(1.55, 1.55, "A", cex=1.2)
draw.circle(1.082293743,1.44944337,0.2,nv=100,border=NULL,col=NA,lty=2,density=NULL,
            angle=45,lwd=1.3)
points(1.082293743,1.44944337,lwd=1.5,col=7)

neighbours = subset(pois2,disc(0.2,c(1.082293743,1.44944337)))
neighbours = subset(neighbours,owin(c(0.25,1.75),c(1.5,2)))
points(neighbours,col=22,main="")

