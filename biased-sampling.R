library(geosphere)
samp <- sample(1:550, size=5)
points <- x[samp,]

distm(x[samp,2:3]/10000, x[samp,2:3]/10000)

apply(distHaversine(points, )
      
apply(c(x, points), 'X', distHaversine)
distHaversine(p1, p2)

dist <- NULL
for(i in 1:dim(points[1])){
  dist[,i] <- distHaversine(points[i,2:3]/10000, x[,2:3]/10000)
}

