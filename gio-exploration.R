setwd('~/Desktop/culex/culex-SDM/')

# load libraries
require(dismo)
require(raster) #needed?
require(maptools)
#require(spatstat)
#require(rworldmap)
#require(scatterplot3d)

# load auxiliary functions
source('code/functions.R')

# baseline environmental data that averages 1950-2000 variables has already been collected in thsi rfile
load('environment-data/africa-baseline/Africa.RData')
resolution <- res(environmental.data)                       # !! we need at least a little info on how the environmental rasters are taken from worldclim data. Can we store only the used parts of the worldclim data on github?

# and so that we can include one IPCC projection:
load('environment-data/africa-forecasts-a2a/Africa-forecast-2050-A2A.RData')
data.a2a <- environmental.data.2050
rm(environmental.data.2050)

ose1 <- rgb(85,108,17, m=255)
ose2 <- rgb(160,108,17, m=255)
ose3 <- rgb(114,132,56, m=255)
ose4 <- rgb(137,152,87, m=255)
# to visualize these colors
palette=c(ose1,ose2,ose3,ose4)
barplot(seq(1:4), col=palette)

mask <- rasterize(wrld_simpl[wrld_simpl$REGION==2,], environmental.data)
environmental.data <- mask(environmental.data, mask)
data.a2a <- mask(data.a2a, mask)
plot(environmental.data[[1]])                                         # !!  what does this plot?
africa.lines <- wrld_simpl[(wrld_simpl$REGION==2 & wrld_simpl$NAME!='South Africa'),]
SA <- wrld_simpl[(wrld_simpl$NAME=='South Africa'),]
SA@polygons[[1]]@Polygons[[2]]@hole<-FALSE # Here we take care of Lesotho  # !! How?
africa.lines2 <- rbind(africa.lines, SA)                        # !! can you just rbind a polygon?
plot(africa.lines2, add=TRUE)

for(i in 1:12) {
name <- paste(month.abb[i], '.temp.diff', sep='')
assign(name, environmental.data[[12+i]]-environmental.data[[i]])
# assign(name, data.a2a[[12+i]]-data.a2a[[i]]) #temperature differential for projections
}

# there are a few other edits to environmental.data, but no need to worry about that now. Notice that the names in the raster brick contain all the things we need

# names(environmental.data)  # shows all of the variables in our set

load('mosq-data/culex.RData')  # this should contain everything I've jsut repeated
par(mar=c(1,1,1,1), mfrow=c(9,10))
for(i in 1:86){
  hist(environmental.data[[i]], xlab='', ylab='', axes=FALSE, main=variable.names[i,1], cex.main=0.7)
  axis(1)
}

pipiens <- read.csv('mosq-data/pipiens.csv', skip=2)
pipiens.coords <- cbind(pipiens$DecimalLongitude, pipiens$DecimalLatitude)
plot(pipiens.coords, pch=20, cex=0.6)

plot(wrld_simpl, add=TRUE, border='grey')  # overlay the world
