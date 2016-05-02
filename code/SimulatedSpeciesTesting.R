library(virtualspecies)
library(dismo)
library(ROCR)
library(spThin)
library(doParallel)
library(foreach)
library(rgeos)
require(sp)
require(spThin)

#157
#123123241
#122345131
N=100

filterByProximity <- function (rec.df.orig, dist, reps=1) #df of x-y, thinning dist, time to repeat thinning process
{
  thin.par <- dist
  reduced.rec.dfs <- list()
  for (Rep in 1:reps) {
    rec.df <- rec.df.orig
    DistMat <- as.matrix(dist(x = rec.df, diag=TRUE))
    diag(DistMat) <- NA #fills in values across the diagonal, if we want to set those to 0 -- use DistMat[upper.tri(DistMat)]=0
    while (min(DistMat, na.rm = TRUE) < thin.par & nrow(rec.df) > 1) {
      CloseRecs <- which(DistMat < thin.par, arr.ind = TRUE)[, 1]
      RemoveRec <- as.numeric(names(which(table(CloseRecs) == max(table(CloseRecs)))))
      if (length(RemoveRec) > 1) {
        RemoveRec <- sample(RemoveRec, 1) #removes clustered points until only one left in the cluster
      }
      rec.df <- rec.df[-RemoveRec, ]
      DistMat <- DistMat[-RemoveRec, -RemoveRec]
      if (length(DistMat) == 1) {
        break
      }
    }
    #colnames(rec.df) <- c("Longitude", "Latitude") #not sure if we want this
    reduced.rec.dfs[[Rep]] <- rec.df
  }
  return(reduced.rec.dfs)
}

time1<-system.time({cl <- makeCluster(2) ;registerDoParallel(cl)

AUC<-foreach (k=1:40, .combine=rbind,.packages=c('virtualspecies', 'dismo','ROCR', 'spThin', 'doParallel', 'rgeos','sp')) %dopar% {
  
set.seed(k)
sp<-generateRandomSp(environmental.data.rs, approach="pca") #generates a species with a random relationship to the two pca axes

set.seed(124)
sample.train<-sampleOccurrences(sp, 2000) # samples 2000 presence points for the species

sample.test<-sampleOccurrences(sp, 2000, type="presence-absence")  ##samples 2000 presence absence points for evaluation

sample.train<-sample.train[[1]][,1:2]
sample.test<-sample.test[[1]][,1:3]

###extract the environmental variables for the training and testing points

train.x<-extract(environmental.data.rs, sample.train)
test.x<-extract(environmental.data.rs, sample.test[,1:2])

###Train models on all of the data

#MAX<-maxent(environmental.data.rs, sample.train)
#Max.full<-predict(MAX, test.x)
#pred<-prediction(Max.full,sample.test[,3])
#AUC.Max<-performance(pred, "auc")@y.values

LOB<-lobag.oc(train.x, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB<-performance(pred, "auc")@y.values

##Establish 3 points (roughly equivalent to the three clusters in our mosquito data) and calculate distances of each training point to each of them
Point1<-c(30,30)
Dist1<-rep(NA, nrow(sample.train))
for (i in 1:nrow(sample.train))
{
  Dist1[i]<-dist(rbind(sample.train[i,], Point1))
}

Point2<-c(10,8)
Dist2<-rep(NA, nrow(sample.train))
for (i in 1:nrow(sample.train))
{
  Dist2[i]<-dist(rbind(sample.train[i,], Point2))
}

Point3<-c(35,5)
Dist3<-rep(NA, nrow(sample.train))
for (i in 1:nrow(sample.train))
{
  Dist3[i]<-dist(rbind(sample.train[i,], Point3))
}

#### Pick the smallest distance to a chosen point for each of the training points
DistM<-cbind(Dist1,Dist2, Dist3)
Distmin<-data.frame(apply(DistM,1,min))
Distmin<-cbind(c(1:nrow(sample.train)), Distmin)

colnames(Distmin)<-c('ID','Distance')

set.seed(345)

### sample from the training points randomly and according to the minimum distance to a chosen point
Rand<-sample(1:nrow(Distmin), size=round(nrow(Distmin))/4, replace=F)
Sample<-sample(1:nrow(Distmin),size=round(nrow(Distmin)/4), replace=F, prob=1/(Distmin$Distance)^30)

RandPoints<-sample.train[Rand,]
SamplePoints<-sample.train[Sample,]

RandX<-train.x[Rand,]
SampleX<-train.x[Sample,]

#plot(sample.train)
#points(RandPoints, col='blue')
#points(SamplePoints, col='red')

################################
####Train models on random, biased, and corrected (corrections below) points and then evaluate at testing points#####
###################################

#MAX<-maxent(environmental.data.rs, RandPoints)
#Max.full<-predict(MAX, test.x)
#pred<-prediction(Max.full,sample.test[,3])
#AUC.Max.Rand<-performance(pred, "auc")@y.values

LOB<-lobag.oc(RandX, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.RAND<-performance(pred, "auc")@y.values

#MAX<-maxent(environmental.data.rs, SamplePoints)
#Max.full<-predict(MAX, test.x)
#pred<-prediction(Max.full,sample.test[,3])
#AUC.Max.Biased<-performance(pred, "auc")@y.values

LOB<-lobag.oc(SampleX, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.Biased<-performance(pred, "auc")@y.values


Grid.points<-gridSample(SamplePoints, environmental.data.rs, n=1)

grid.x<-extract(environmental.data.rs, Grid.points)


#MAX<-maxent(environmental.data.rs, Grid.points)
#Max.full<-predict(MAX, test.x)
#pred<-prediction(Max.full,sample.test[,3])
#AUC.Max.Grid<-performance(pred, "auc")@y.values

LOB<-lobag.oc(grid.x, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.Grid<-performance(pred, "auc")@y.values

res=1
r <- raster(extent(range(SamplePoints[,1]), range(SamplePoints[,2])) + res)
res(r) <- res

Grid.points1<-gridSample(SamplePoints, r, n=1)

grid.x1<-extract(environmental.data.rs, Grid.points1)

#MAX<-maxent(environmental.data.rs, Grid.points1)
#Max.full<-predict(MAX, test.x)
#pred<-prediction(Max.full,sample.test[,3])
#AUC.Max.Grid1<-performance(pred, "auc")@y.values

LOB<-lobag.oc(grid.x1, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.Grid1<-performance(pred, "auc")@y.values

res=5
r <- raster(extent(range(SamplePoints[,1]), range(SamplePoints[,2])) + res)
res(r) <- res

Grid.points5<-gridSample(SamplePoints, r, n=1)


grid.x5<-extract(environmental.data.rs, Grid.points5)


LOB<-lobag.oc(grid.x5, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.Grid5<-performance(pred, "auc")@y.values


####Distance spatial thinning##############

pts <- as.matrix(SamplePoints)
pts2 <- filterByProximity(pts,dist=.167)
pts.ind<-SampleX[which(pts[,1] %in% pts2[[1]][,1]),]

LOB<-lobag.oc(pts.ind, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.dist.16<-performance(pred, "auc")@y.values

pts <- as.matrix(SamplePoints)
pts2 <- filterByProximity(pts,dist=1)
pts.ind<-SampleX[which(pts[,1] %in% pts2[[1]][,1]),]

LOB<-lobag.oc(pts.ind, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.dist1<-performance(pred, "auc")@y.values

pts <- as.matrix(SamplePoints)
pts2 <- filterByProximity(pts,dist=5)
pts.ind<-SampleX[which(pts[,1] %in% pts2[[1]][,1]),]

LOB<-lobag.oc(pts.ind, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.dist5<-performance(pred, "auc")@y.values

# Point1<-c(30,30)
# Dist1<-rep(NA, nrow(SamplePoints))
# for (i in 1:nrow(SamplePoints))
# {
#   Dist1[i]<-dist(rbind(SamplePoints[i,], Point1))
# }
# 
# 
# Point2<-c(10,8)
# Dist2<-rep(NA, nrow(SamplePoints))
# for (i in 1:nrow(SamplePoints))
# {
#   Dist2[i]<-dist(rbind(SamplePoints[i,], Point2))
# }
# 
# Point3<-c(35,5)
# Dist3<-rep(NA, nrow(SamplePoints))
# for (i in 1:nrow(SamplePoints))
# {
#   Dist3[i]<-dist(rbind(SamplePoints[i,], Point3))
# }
# 
# 
# 
# #### Pick the smallest distance to a chosen point for each of the training points
# DistM<-cbind(Dist1,Dist2, Dist3)
# Distmin<-data.frame(apply(DistM,1,min))
# Distmin<-cbind(c(1:nrow(SamplePoints)), Distmin)
# 
# 
# colnames(Distmin)<-c('ID','Distance')
# 
# 
# 
# resampled<-sample(1:nrow(SampleX), size=nrow(SampleX), replace=T, prob=(Distmin[,2]^30))
# resampled<-SampleX[resampled,]
# 
# 
# LOB<-lobag.oc(resampled, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.RESAMP<-performance(pred, "auc")@y.values
# 
# resampled<-sample(1:nrow(SampleX), size=2*nrow(SampleX)/3, replace=T, prob=(Distmin[,2]^30))
# resampled<-SampleX[resampled,]
# 
# 
# LOB<-lobag.oc(resampled, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.RESAMP.66<-performance(pred, "auc")@y.values
# 
# resampled<-sample(1:nrow(SampleX), size=nrow(SampleX)/2, replace=T, prob=(Distmin[,2]^30))
# resampled<-SampleX[resampled,]
# 
# 
# LOB<-lobag.oc(resampled, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.RESAMP.5<-performance(pred, "auc")@y.values
# 
# 
# resampled<-sample(1:nrow(SampleX), size=nrow(SampleX)/3, replace=T, prob=(Distmin[,2]^30))
# resampled<-SampleX[resampled,]
# 
# 
# LOB<-lobag.oc(resampled, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.RESAMP.33<-performance(pred, "auc")@y.values
# 
# resampled<-sample(1:nrow(SampleX), size=nrow(SampleX)/4, replace=T, prob=(Distmin[,2]^30))
# resampled<-SampleX[resampled,]
# 
# 
# LOB<-lobag.oc(resampled, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.RESAMP.25<-performance(pred, "auc")@y.values

###################################Environmental Thinning#################################


unique.pca<-predict(pca,SampleX)
ID <- c(1:length(unique.pca[,1]))
unique.pca <- cbind(ID, unique.pca)
plot(unique.pca[,2:3])

# #assign raster layer on top of unique pca to sample from
# r<-raster(xmn=min(unique.pca[,2]-1), xmx=max(unique.pca[,2]+1), ymn=min(unique.pca[,3]-1), ymx=max(unique.pca[,3]+1))
# res(r)<-(.25) #size of grid squares
# cell<-cellFromXY(r, unique.pca[,2:3])
# dup<-duplicated(cell)
# data.thin<- unique.pca[!dup,]
# plot(data.thin[,2:3])
# #coordinate in geographic space
# pca.coords <- SamplePoints[data.thin[,1],] #coordinates of pca points 
# pca.env<-extract(environmental.data.rs, pca.coords) #how is this right?? len(afr.pip.unique)=430 but unique.x was only 402, so how do they match up?
# 
# LOB<-lobag.oc(pca.env, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.pca.grid.25<-performance(pred, "auc")@y.values

r<-raster(xmn=min(unique.pca[,2]-1), xmx=max(unique.pca[,2]+1), ymn=min(unique.pca[,3]-1), ymx=max(unique.pca[,3]+1))
res(r)<-(.5) #size of grid squares
cell<-cellFromXY(r, unique.pca[,2:3])
dup<-duplicated(cell)
data.thin<- unique.pca[!dup,]
plot(data.thin[,2:3])
#coordinate in geographic space
pca.coords <- SamplePoints[data.thin[,1],] #coordinates of pca points 
pca.env<-extract(environmental.data.rs, pca.coords) #how is this right?? len(afr.pip.unique)=430 but unique.x was only 402, so how do they match up?

LOB<-lobag.oc(pca.env, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.pca.grid.5<-performance(pred, "auc")@y.values

# r<-raster(xmn=min(unique.pca[,2]-1), xmx=max(unique.pca[,2]+1), ymn=min(unique.pca[,3]-1), ymx=max(unique.pca[,3]+1))
# res(r)<-(.75) #size of grid squares
# cell<-cellFromXY(r, unique.pca[,2:3])
# dup<-duplicated(cell)
# data.thin<- unique.pca[!dup,]
# plot(data.thin[,2:3])
# #coordinate in geographic space
# pca.coords <- SamplePoints[data.thin[,1],] #coordinates of pca points 
# pca.env<-extract(environmental.data.rs, pca.coords) #how is this right?? len(afr.pip.unique)=430 but unique.x was only 402, so how do they match up?
# 
# LOB<-lobag.oc(pca.env, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.pca.grid.75<-performance(pred, "auc")@y.values


r<-raster(xmn=min(unique.pca[,2]-1), xmx=max(unique.pca[,2]+1), ymn=min(unique.pca[,3]-1), ymx=max(unique.pca[,3]+1))
res(r)<-(1) #size of grid squares
cell<-cellFromXY(r, unique.pca[,2:3])
dup<-duplicated(cell)
data.thin<- unique.pca[!dup,]
plot(data.thin[,2:3])
#coordinate in geographic space
pca.coords <- SamplePoints[data.thin[,1],] #coordinates of pca points 
pca.env<-extract(environmental.data.rs, pca.coords) #how is this right?? len(afr.pip.unique)=430 but unique.x was only 402, so how do they match up?


LOB<-lobag.oc(pca.env, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.pca.grid1<-performance(pred, "auc")@y.values

r<-raster(xmn=min(unique.pca[,2]-1), xmx=max(unique.pca[,2]+1), ymn=min(unique.pca[,3]-1), ymx=max(unique.pca[,3]+1))
res(r)<-(2) #size of grid squares
cell<-cellFromXY(r, unique.pca[,2:3])
dup<-duplicated(cell)
data.thin<- unique.pca[!dup,]
plot(data.thin[,2:3])
#coordinate in geographic space
pca.coords <- SamplePoints[data.thin[,1],] #coordinates of pca points 
pca.env<-extract(environmental.data.rs, pca.coords) #how is this right?? len(afr.pip.unique)=430 but unique.x was only 402, so how do they match up?

LOB<-lobag.oc(pca.env, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.pca.grid2<-performance(pred, "auc")@y.values


#############Distance environmental thinning################


#function for proximity filtering

# pts <- unique.pca[,2:3]
# pts2 <- filterByProximity(pts,dist=.25)
# pts.ind<-SampleX[which(pts[,1] %in% pts2[[1]][,1]),]
# 
# 
# LOB<-lobag.oc(pts.ind, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.pca.dist.25<-performance(pred, "auc")@y.values

pts <- unique.pca[,2:3]
pts2 <- filterByProximity(pts,dist=.5)
pts.ind<-SampleX[which(pts[,1] %in% pts2[[1]][,1]),]


LOB<-lobag.oc(pts.ind, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.pca.dist.5<-performance(pred, "auc")@y.values

# pts <- unique.pca[,2:3]
# pts2 <- filterByProximity(pts,dist=.75)
# pts.ind<-SampleX[which(pts[,1] %in% pts2[[1]][,1]),]
# 
# 
# LOB<-lobag.oc(pts.ind, n.votes=250)
# LOB.full<-predictSvm(LOB,test.x)
# pred<-prediction(LOB.full$p.out,sample.test[,3])
# AUC.LOB.pca.dist.75<-performance(pred, "auc")@y.values

pts <- unique.pca[,2:3]
pts2 <- filterByProximity(pts,dist=1)
pts.ind<-SampleX[which(pts[,1] %in% pts2[[1]][,1]),]


LOB<-lobag.oc(pts.ind, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.pca.dist1<-performance(pred, "auc")@y.values

pts <- unique.pca[,2:3]
pts2 <- filterByProximity(pts,dist=2)
pts.ind<-SampleX[which(pts[,1] %in% pts2[[1]][,1]),]


LOB<-lobag.oc(pts.ind, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.pca.dist2<-performance(pred, "auc")@y.values



AUCs<-cbind(as.numeric(AUC.LOB),as.numeric(AUC.LOB.RAND),as.numeric(AUC.LOB.Biased),
            as.numeric(AUC.LOB.Grid),as.numeric(AUC.LOB.Grid1),as.numeric(AUC.LOB.Grid5), 
            as.numeric(AUC.LOB.dist.16), as.numeric(AUC.LOB.dist1),as.numeric(AUC.LOB.dist5),
            as.numeric(AUC.LOB.pca.grid.5),
            as.numeric(AUC.LOB.pca.grid1),as.numeric(AUC.LOB.pca.grid2),
            as.numeric(AUC.LOB.pca.dist.5),
            as.numeric(AUC.LOB.pca.dist1),as.numeric(AUC.LOB.pca.dist2))
return(AUCs)
}})
stopCluster(cl)
colnames(AUC)<-c('AUC.LOB','AUC.LOB.RAND','AUC.LOB.Biased','AUC.LOB.Grid','AUC.LOB.Grid1','AUC.LOB.Grid5','AUC.LOB.dist.16','AUC.LOB.dist1','AUC.LOB.dist5', 'AUC.LOB.pca.grid.5','AUC.LOB.pca.grid1','AUC.LOB.pca.grid2','AUC.LOB.pca.dist.5','AUC.LOB.pca.dist1','AUC.LOB.pca.dist2')

boxplot(AUC)

Dfunc<-function(x) (x[4:15]-x[3])/(x[1]-x[3])
 DAUC<-t(apply(AUC, MARGIN=1, FUN=Dfunc))
 
 BestModel<-function(x) which(x==max(x))
 
 Best<-apply(DAUC, MARGIN=1, FUN=BestModel)
 Best<-colnames(DAUC)[Best]
 
 DAUCvec<-data.frame(cbind(as.numeric(as.vector(DAUC)),c(rep("Geographic", 180), rep("Environmental", 180)),c(rep("Grid", 90), rep("Distance", 90), rep("Grid", 90), rep("Distance", 90)),c(rep(".16", 30),rep('1',30),rep('5',30),rep(".16", 30),rep('1',30),rep('5',30), rep('.5',30), rep('1', 30), rep('2', 30),rep('.5',30), rep('1', 30), rep('2', 30)),c(rep(colnames(DAUC)[1],30),rep(colnames(DAUC)[2],30),rep(colnames(DAUC)[3],30),rep(colnames(DAUC)[4],30),rep(colnames(DAUC)[5],30),rep(colnames(DAUC)[6],30),rep(colnames(DAUC)[7],30),rep(colnames(DAUC)[8],30),rep(colnames(DAUC)[9],30),rep(colnames(DAUC)[10],30),rep(colnames(DAUC)[11],30),rep(colnames(DAUC)[12],30)), rep(seq(1,30),12)))
 
 
 colnames(DAUCvec)<-c("DAUC", "Space","Type","Distance","Method", "Species")
 
 ModelFull<-lm(as.numeric(DAUC)~as.factor(Space) + as.factor(Type) + as.factor(Distance), data=DAUCvec)
 summary(ModelFull)
 summary(aov(ModelFull))
 
 tukey<-TukeyHSD(aov(ModelFull))
 tukey$`DAUCvec[, 2]`[which(tukey$`DAUCvec[, 2]`[,4]<.05),1]

 summary(aov(as.numeric(DAUCvec[,1])~DAUCvec[,2]+DAUCvec[,3]))
 tukey<-TukeyHSD(aov(as.numeric(DAUCvec[,1])~DAUCvec[,2]))
 tukey$`DAUCvec[, 2]`[which(tukey$`DAUCvec[, 2]`[,4]<.05),1]
 summary(lm(as.numeric(DAUCvec[,1])~DAUCvec[,2]+DAUCvec[,3]))
 
 levelplot(DAUC, xlab="Species")

meanAUC<-apply(AUC, MARGIN=2, FUN=mean)

#AUCs<-data.frame(cbind(as.numeric(AUC.LOB),as.numeric(AUC.LOB.RAND),as.numeric(AUC.LOB.Biased),as.numeric(AUC.LOB.Grid),as.numeric(AUC.LOB.Grid1),as.numeric(AUC.LOB.Grid5), as.numeric(AUC.LOB.pca.grid.25),as.numeric(AUC.LOB.pca.grid.5),as.numeric(AUC.LOB.pca.grid.75),as.numeric(AUC.LOB.pca.grid1),as.numeric(AUC.LOB.pca.grid2)))

#as.numeric(AUC.LOB.RESAMP),as.numeric(AUC.LOB.RESAMP.66),as.numeric(AUC.LOB.RESAMP.5), as.numeric(AUC.LOB.RESAMP.33), as.numeric(AUC.LOB.RESAMP.25)

###Find mean AUCs across all random species
#meanAUCs<-apply(AUCs, MARGIN=2, FUN=mean)

###calculate delta AUC according to Fourcade et al. 
DAUCmean<-colMeans(DAUC)

DAUCmean2<-(meanAUC[4:15]-meanAUC[3])/(meanAUC[1]-meanAUC[3])
boxplot(DAUC)

save.image("SimulatedSpeciesTesting")
