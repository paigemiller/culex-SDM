library(dismo)
library(kernlab)
### data with PCA data included
load("/home/pbmpb13/Documents/culex-SDM/mosq-data/culex-v4.RData")
afr.pip<-pipiens.coords[which(pipiens.coords[,1]>-35&pipiens.coords[,1]<38&pipiens.coords[,2]>-20&pipiens.coords[,2]<50),]
plot(wrld_simpl)
points(afr.pip, col='red')

###subsetting data by unique coordinates
afr.pip.unique<-unique(afr.pip)
plot(africa.lines2)
points(afr.pip.unique, col="red")

#data frame of xy points
afr.pip.unique<-data.frame(afr.pip.unique)

#environmental data at unique points in geographical space
unique.x<-extract(environmental.data.rs, afr.pip.unique)
unique.x<-data.frame(unique.x[complete.cases(unique.x),]) #take out NAs


###Basic model removing duplicate points
LOB<-lobag.oc(unique.x, n.votes=250)
#sum votes
rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  rastersum<-rastersum+rasterpredict
}
rastermean.unique<-rastersum/length(LOB)

plot(rastermean.unique)

### Environmental grid thin
unique.pca<-predict(pca,unique.x)
ID <- c(1:length(unique.pca[,1]))
unique.pca <- cbind(ID, unique.pca)
plot(unique.pca[,2:3])

#assign raster layer on top of unique pca to sample from 
r<-raster(xmn=min(unique.pca[,2]-1), xmx=max(unique.pca[,2]+1), ymn=min(unique.pca[,3]-1), ymx=max(unique.pca[,3]+1))
res(r)<-(.1) #size of grid squares
cell<-cellFromXY(r, unique.pca[,2:3])
dup<-duplicated(cell)
data.thin<- unique.pca[!dup,]
plot(data.thin[,2:3])
#coordinate in geographic space
pca.coords <- afr.pip.unique[data.thin[,1],] #coordinates of pca points 
pca.env<-extract(environmental.data.rs, pca.coords) #how is this right?? len(afr.pip.unique)=430 but unique.x was only 402, so how do they match up?

#Lobag model
LOB.pca.env<-lobag.oc(pca.env, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB.pca.env[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB.pca.env[[i]])
  rastersum<-rastersum+rasterpredict
}
rastermean.pca.env<-rastersum/length(LOB)

plot(rastermean.pca.env)


### Environmental distance thinning
#pca.dist <- pointDistance(unique.pca[,2:3], lonlat=F)

#function for proximity filtering
filterByProximity <- function(xy, dist, mapUnits = F) {
  #xy can be either a SpatialPoints or SPDF object, or a matrix
  #dist is in km if mapUnits=F, in mapUnits otherwise
  if (!mapUnits) {
    d <- spDists(xy,longlat=T)
  }
  if (mapUnits) {
    d <- spDists(xy[,2:3],longlat=F)
  }
  diag(d) <- NA
  close <- (d <= dist)
  diag(close) <- NA
  closePts <- which(close,arr.ind=T)
  discard <- matrix(nrow=2,ncol=2)
  if (nrow(closePts) > 0) {
    while (nrow(closePts) > 0) {
      if ((!paste(closePts[1,1],closePts[1,2],sep='_') %in% paste(discard[,1],discard[,2],sep='_')) & (!paste(closePts[1,2],closePts[1,1],sep='_') %in% paste(discard[,1],discard[,2],sep='_'))) {
        discard <- rbind(discard, closePts[1,])
        closePts <- closePts[-union(which(closePts[,1] == closePts[1,1]), which(closePts[,2] == closePts[1,1])),]
      }
    }
    discard <- discard[complete.cases(discard),]
    return(xy[-discard[,1],])
  }
  if (nrow(closePts) == 0) {
    return(xy)
  }
}

library(rgeos)
require(sp)
pts <- unique.pca[,2:3]
pts2 <- filterByProximity(pts,dist=2, mapUnits=T)

plot(pts)
axis(1)
axis(2)
apply(as.data.frame(pts),1,function(x) plot(gBuffer(SpatialPoints(coords=matrix(c(x[1],x[2]),nrow=1)),width=2),add=T))
par(new=T)
plot(pts2,add=T,col='blue',pch=20,cex=2)

pts.ind<-unique.x[which(pts2[,1] %in% pts[,1]),]

#Lobag model for environmental distance thin
LOB.pca.dist<-lobag.oc(pts.ind, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB.pca.dist[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB.pca.dist[[i]])
  rastersum<-rastersum+rasterpredict
}
rastermean.pca.dist<-rastersum/length(LOB)

plot(rastermean.pca.dist)

