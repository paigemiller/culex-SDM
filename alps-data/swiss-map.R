library("ggmap")
library(maptools)
library(maps)
library(mapproj)

## environmental sampling
dat <- data.frame (lat=x$X/10000, lon=x$Y/10000)

map <- get_map(location = 'Switzerland', zoom = 7, maptype="roadmap")
ggmap(map) + 
  geom_point(aes(x=lat, y=lon), alpha = .5, color="darkred")

