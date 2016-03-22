load("C:/Users/Robert/OneDrive/Documents/Robbie/8910ENM/mosq-data/culex-v4.RData")
library(spThin)
library(dismo)
afr.pip<-pipiens.coords[which(pipiens.coords[,1]>-35&pipiens.coords[,1]<38&pipiens.coords[,2]>-20&pipiens.coords[,2]<50),]
plot(wrld_simpl)
points(afr.pip, col='red')

afr.pip.unique<-unique(afr.pip)
plot(africa.lines)
points(afr.pip.unique)

afr.pip.unique<-data.frame(afr.pip.unique)
SPEC<-rep(1,nrow(afr.pip.unique))
thinned<-thin(cbind(afr.pip.unique,SPEC), lat.col="X1", long.col="X2", spec.col='SPEC', thin.par=10, reps=100, locs.thinned.list.return = T, write.files = F, verbose=T)

plot(africa.lines)
plot(wrld_simpl)
points(thinned[[3]][,2],thinned[[3]][,1], col="red")
nrow(thinned[[3]])
thinned[[2]]
thinned.points<-cbind(thinned[[3]][,2], thinned[[3]][,1])
thinned.x<-extract(environmental.data.rs, thinned.points)
thinned.x<-thinned.x[complete.cases(thinned.x),]

unique.x<-extract(environmental.data.rs, afr.pip.unique)
unique.x<-unique.x[complete.cases(unique.x),]


LOB<-lobag.oc(thinned.x, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.thinned<-rastersum/length(LOB)

plot(rastermean.thinned)

LOB<-lobag.oc(unique.x, n.votes=250)

rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.unique<-rastersum/length(LOB)

plot(rastermean.unique)


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
