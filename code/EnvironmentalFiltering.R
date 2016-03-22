load("C:/Users/Robert/OneDrive/Documents/Robbie/8910ENM/mosq-data/culex-v4.RData")
library(dismo)
afr.pip<-pipiens.coords[which(pipiens.coords[,1]>-35&pipiens.coords[,1]<38&pipiens.coords[,2]>-20&pipiens.coords[,2]<50),]
#plot(wrld_simpl)
#points(afr.pip, col='red')

afr.pip.unique<-unique(afr.pip)
plot(africa.lines2)
points(afr.pip.unique)

afr.pip.unique<-data.frame(afr.pip.unique)



unique.x<-extract(environmental.data.rs, afr.pip.unique)
unique.x<-unique.x[complete.cases(unique.x),]

unique.pca<-predict(pca,unique.x)

plot(unique.pca)
r<-raster(xmn=min(unique.pca[,1]-1), xmx=max(unique.pca[,1]+1), ymn=min(unique.pca[,2]-1), ymx=max(unique.pca[,2]+1))
res(r)<-(1)
cell<-cellFromXY(r, unique.pca)
dup<-duplicated(cell)
data.thin<- unique.pca[!dup,]
plot(data.thin[,1:2])


###Unique points
LOB<-lobag.oc(unique.x, n.votes=250)

rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.unique<-rastersum/length(LOB)

plot(rastermean.unique)


## Distance Thinned
LOB<-lobag.oc(thinned.x, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.thinned<-rastersum/length(LOB)

plot(rastermean.thinned)


### Grid Thinned
Grid.points<-gridSample(afr.pip.unique, environmental.data.rs, n=1)

grid.x<-extract(environmental.data.rs, Grid.points)
grid.x<-grid.x[complete.cases(grid.x),]

LOB<-lobag.oc(grid.x, n.votes=250)

rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.grid<-rastersum/length(LOB)

plot(rastermean.grid)
