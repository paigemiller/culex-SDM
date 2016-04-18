## script for culex quinq

afr.quq<-quinq.coords[which(quinq.coords[,1]>-35&quinq.coords[,1]<38&quinq.coords[,2]>-20&quinq.coords[,2]<50),]
plot(wrld_simpl)
points(afr.quq, col='red')

afr.quq.unique<-unique(afr.quq)
plot(africa.lines)
points(afr.quq.unique, col="red")

afr.pip.unique<-data.frame(afr.pip.unique)

unique.x<-extract(environmental.data.rs, afr.quq.unique)
unique.x<-unique.x[complete.cases(unique.x),]

SPEC<-rep(1,nrow(afr.quq.unique))
thinned<-thin(cbind(afr.quq.unique,SPEC), lat.col="X1", long.col="X2", spec.col='SPEC', thin.par=10, 
              reps=100, locs.thinned.list.return = T, write.files = F, verbose=T)

thinned100<-thin(cbind(afr.pip.unique,SPEC), lat.col="X1", long.col="X2", spec.col='SPEC', thin.par=100, 
                 reps=100, locs.thinned.list.return = T, write.files = F, verbose=T)

thinned500<-thin(cbind(afr.quq.unique,SPEC), lat.col="X1", long.col="X2", spec.col='SPEC', thin.par=500, 
                 reps=100, locs.thinned.list.return = T, write.files = F, verbose=T)

pip.pop<-extract(pop, afr.pip.unique)

pip.pop<-(pip.pop+1)/max(pip.pop)

pip.resample<-sample(x=1:nrow(afr.pip.unique),size=nrow(afr.pip.unique),replace=T, prob=(1/(pip.pop)))

resample.points<-afr.pip.unique[pip.resample,]

resample.x<-extract(environmental.data.rs, resample.points)
resample.x<-resample.x[complete.cases(resample.x),]

#plot(africa.lines)
#plot(wrld_simpl)
#points(thinned[[3]][,2],thinned[[3]][,1], col="red")
nrow(thinned[[5]])
#thinned[[2]]
thinned.points<-cbind(thinned[[5]][,2], thinned[[5]][,1])
thinned.x<-extract(environmental.data.rs, thinned.points)
thinned.x<-thinned.x[complete.cases(thinned.x),]

thinned.points100<-cbind(thinned100[[5]][,2], thinned100[[5]][,1])
thinned.x100<-extract(environmental.data.rs, thinned.points100)
thinned.x100<-thinned.x100[complete.cases(thinned.x100),]

thinned.points500<-cbind(thinned500[[5]][,2], thinned500[[5]][,1])
thinned.x500<-extract(environmental.data.rs, thinned.points500)
thinned.x500<-thinned.x500[complete.cases(thinned.x500),]

Grid.points<-gridSample(afr.pip.unique, environmental.data.rs, n=1)

grid.x<-extract(environmental.data.rs, Grid.points)
grid.x<-grid.x[complete.cases(grid.x),]

res=1
r <- raster(extent(range(afr.pip.unique[,1]), range(afr.pip.unique[,2])) + res)
res(r) <- res

Grid.points1<-gridSample(afr.pip.unique, r, n=1)

grid.x1<-extract(environmental.data.rs, Grid.points1)
grid.x1<-grid.x1[complete.cases(grid.x1),]

res=5
r <- raster(extent(range(afr.pip.unique[,1]), range(afr.pip.unique[,2])) + res)
res(r) <- res

Grid.points5<-gridSample(afr.pip.unique, r, n=1)

grid.x5<-extract(environmental.data.rs, Grid.points5)
grid.x5<-grid.x5[complete.cases(grid.x5),]

###Environmental Grid
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

pts <- unique.pca[,2:3]
pts2 <- filterByProximity(pts,dist=.1, mapUnits=F)

plot(pts)
axis(1)
axis(2)
#apply(as.data.frame(pts),1,function(x) plot(gBuffer(SpatialPoints(coords=matrix(c(x[1],x[2]),nrow=1)),width=2),add=T))
par(new=T)
plot(pts2,add=T,col='blue',pch=20,cex=2)

pts.ind<-unique.x[which(pts2[,1] %in% pts[,1]),]
```

```{r, eval=F, echo=F}
MAX<-maxent(environmental.data.rs, thinned.points100)
Max.thinned100<-predict(environmental.data.rs, MAX)

plot(Max.thinned100)


LOB<-lobag.oc(thinned.x100, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.thinned100<-rastersum/length(LOB)

plot(rastermean.thinned100)


MAX<-maxent(environmental.data.rs, thinned.points500)
Max.thinned500<-predict(environmental.data.rs, MAX)

plot(Max.thinned500)

points(thinned.points500)

LOB<-lobag.oc(thinned.x500, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.thinned500<-rastersum/length(LOB)

plot(rastermean.thinned500)
points(thinned.points500)

MAX<-maxent(environmental.data.rs, Grid.points1)
Max.grid1<-predict(environmental.data.rs, MAX)

plot(Max.grid1)

LOB<-lobag.oc(grid.x1, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.grid1<-rastersum/length(LOB)

plot(rastermean.grid1)


MAX<-maxent(environmental.data.rs, Grid.points5)
Max.grid5<-predict(environmental.data.rs, MAX)

plot(Max.grid5)


LOB<-lobag.oc(grid.x5, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.grid5<-rastersum/length(LOB)

plot(rastermean.grid5)
```


```{r LoBag Models, echo=F, eval=F}
LOB<-lobag.oc(unique.x, n.votes=250)

rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.unique<-rastersum/length(LOB)

plot(rastermean.unique)

MAX<-maxent(environmental.data.rs, afr.pip.unique)
Max.unique<-predict(environmental.data.rs, MAX)

##################################################
LOB<-lobag.oc(thinned.x, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.thinned<-rastersum/length(LOB)

plot(rastermean.thinned)

MAX<-maxent(environmental.data.rs, thinned.points)
Max.thinned<-predict(environmental.data.rs,MAX)


##############################################
LOB<-lobag.oc(grid.x, n.votes=250)

rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.grid<-rastersum/length(LOB)

MAX<-maxent(environmental.data.rs,Grid.points)
Max.grid<-predict(environmental.data.rs, MAX)

plot(rastermean.grid)

###########################################

LOB<-lobag.oc(resample.x, n.votes=250)

rastersum<- predict(environmental.data.rs,LOB[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB[[i]])
  
  rastersum<-rastersum+rasterpredict
}
rastermean.resample<-rastersum/length(LOB)

MAX<-maxent(environmental.data.rs, resample.points)
Max.resample<-predict(environmental.data.rs,MAX)

plot(rastermean.resample)

#Grid Thinned Environmental
LOB.pca.env<-lobag.oc(pca.env, n.votes=250)
rastersum<- predict(environmental.data.rs,LOB.pca.env[[1]])
for (i in 2:250)
{
  rasterpredict<-predict(environmental.data.rs,LOB.pca.env[[i]])
  rastersum<-rastersum+rasterpredict
}
rastermean.pca.env<-rastersum/length(LOB)

plot(rastermean.pca.env)

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


save.image('Spatial_Thinning')

```