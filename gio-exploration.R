setwd('~/Desktop/culex/culex-SDM/')

# load libraries
require(dismo)
require(raster) 
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

load('mosq-data/culex-v2.RData')  # this should contain everything I've just repeated, but there are four versions
par(mar=c(1,1,1,1), mfrow=c(9,10))
for(i in 1:86){
  hist(environmental.data[[i]], xlab='', ylab='', axes=FALSE, main=variable.names[i,1], cex.main=0.7)
  axis(1)
}

# load Cx. pipiens data
pipiens <- read.csv('mosq-data/pipiens.csv', skip=2)
pipiens.coords <- cbind(pipiens$DecimalLongitude, pipiens$DecimalLatitude)

quinq <- read.csv('mosq-data/quinq.csv', skip=2)
quinq.coords <- cbind(quinq$DecimalLongitude, quinq$DecimalLatitude)

salinarius <- read.csv('mosq-data/salinarius.csv', skip=2)
salinarius.coords <- cbind(salinarius$DecimalLongitude, salinarius$DecimalLatitude)

plot(quinq.coords, pch=20, cex=0.6, xlim=c(-25,55))
points(pipiens.coords, xlim=c(-25,55), pch=21, cex=0.6)
plot(wrld_simpl, add=TRUE, border='grey')  # overlay the world

# overlaying pipiens presence points on January min temp
plot(environmental.data.rs[[1]])
points(pipiens.coords, pch=21, bg='yellow', cex=0.6)
# getting the data for just those points
pipiens.presence <- data.frame(extract(environmental.data.rs, pipiens.coords))
z <- which(is.na(pipiens.presence[,1]))
pipiens.presence <- pipiens.presence[-z,]  # now the env data for cells where coords are found
n.presence <- dim(pipiens.presence)[1]
pts.presence <- pipiens.coords[-z,]  # and corresponding coordinates

# simulate background data
set.seed(3281994)
backgr.pts <- randomPoints(environmental.data.rs, n.presence)
data.background <- data.frame(extract(environmental.data.rs, backgr.pts))  # !! does this give the 86-wide vector of environmental data for each raster cell that backgr.pts is in?
cols <- colorRampPalette(c('darkorange4', 'khaki1'), interpolate='linear')

# plotting points on background of bioclim
plot(environmental.data[[37]]/10, col=(cols(100)), main='Presence/background points')
points(pipiens.coords, pch=21, bg='yellow', cex=0.6)
n <- sample(n.presence, 100)
points(backgr.pts[n,], pch=21, bg='grey', cex=0.6) #  too many to plot
plot(wrld_simpl[wrld_simpl$REGION==2,], add=TRUE)
legend('topright', pch=21, pt.bg=c('yellow','grey'), legend=c('Cx. pipiens', 'Background'), bty='n', cex=0.7)
text(-17,-25,'Annual mean', pos=4, cex=0.8)
text(-17,-28, 'temperature (degrees C)', pos=4, cex=0.8)

# culex-v3.RData contains up to this point

# split training and testing datasets 
set.seed(3281994)
train.presence.id <- sample(seq(1,n.presence),ceiling(0.8*n.presence))
train.background.id <- sample(seq(1,n.presence),ceiling(0.8*n.presence))
train.presence <- pipiens.presence[train.presence.id,]
train.background <- data.background[train.background.id,]
test.presence <- pipiens.presence[-train.presence.id,]
test.background <- data.background[-train.presence.id,]
