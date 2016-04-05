library(virtualspecies)
library(dismo)
library(ROCR)
library(spThin)
#157
#123123241
#122345131



AUC.LOB<-rep(NA,30)
AUC.LOB.RAND<-rep(NA,30)
AUC.LOB.Biased<-rep(NA,30)
AUC.LOB.Grid<-rep(NA,30)
AUC.LOB.Grid1<-rep(NA,30)
AUC.LOB.Grid5<-rep(NA,30)
AUC.LOB.pca.grid<-rep(NA,30)

for (k in 1:30)
{
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
AUC.LOB[k]<-performance(pred, "auc")@y.values



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
AUC.LOB.RAND[k]<-performance(pred, "auc")@y.values



#MAX<-maxent(environmental.data.rs, SamplePoints)
#Max.full<-predict(MAX, test.x)
#pred<-prediction(Max.full,sample.test[,3])
#AUC.Max.Biased<-performance(pred, "auc")@y.values

LOB<-lobag.oc(SampleX, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.Biased[k]<-performance(pred, "auc")@y.values


Grid.points<-gridSample(SamplePoints, environmental.data.rs, n=1)

grid.x<-extract(environmental.data.rs, Grid.points)


#MAX<-maxent(environmental.data.rs, Grid.points)
#Max.full<-predict(MAX, test.x)
#pred<-prediction(Max.full,sample.test[,3])
#AUC.Max.Grid<-performance(pred, "auc")@y.values

LOB<-lobag.oc(grid.x, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.Grid[k]<-performance(pred, "auc")@y.values

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
AUC.LOB.Grid1[k]<-performance(pred, "auc")@y.values

res=5
r <- raster(extent(range(SamplePoints[,1]), range(SamplePoints[,2])) + res)
res(r) <- res

Grid.points5<-gridSample(SamplePoints, r, n=1)


grid.x5<-extract(environmental.data.rs, Grid.points5)


LOB<-lobag.oc(grid.x5, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.Grid5[k]<-performance(pred, "auc")@y.values






###################################Environmental Thinning#################################


unique.pca<-predict(pca,SampleX)
ID <- c(1:length(unique.pca[,1]))
unique.pca <- cbind(ID, unique.pca)
plot(unique.pca[,2:3])

#assign raster layer on top of unique pca to sample from 
r<-raster(xmn=min(unique.pca[,2]-1), xmx=max(unique.pca[,2]+1), ymn=min(unique.pca[,3]-1), ymx=max(unique.pca[,3]+1))
res(r)<-(1) #size of grid squares
cell<-cellFromXY(r, unique.pca[,2:3])
dup<-duplicated(cell)
data.thin<- unique.pca[!dup,]
plot(data.thin[,2:3])
#coordinate in geographic space
pca.coords <- SamplePoints[data.thin[,1],] #coordinates of pca points 
pca.env<-extract(environmental.data.rs, pca.coords) #how is this right?? len(afr.pip.unique)=430 but unique.x was only 402, so how do they match up?





#############################




LOB<-lobag.oc(pca.env, n.votes=250)
LOB.full<-predictSvm(LOB,test.x)
pred<-prediction(LOB.full$p.out,sample.test[,3])
AUC.LOB.pca.grid[k]<-performance(pred, "auc")@y.values

}



AUCs<-data.frame(cbind(as.numeric(AUC.LOB),as.numeric(AUC.LOB.RAND),as.numeric(AUC.LOB.Biased),as.numeric(AUC.LOB.Grid),as.numeric(AUC.LOB.Grid1),as.numeric(AUC.LOB.Grid5), as.numeric(AUC.LOB.pca.grid)))


###Find mean AUCs across all random species
meanAUCs<-apply(AUCs, MARGIN=2, FUN=mean)

###calculate delta AUC according to Fourcade et al. 
DAUC<-meanAUCs[4:7]

DAUC<-(DAUC-meanAUCs[3])/(meanAUCs[1]-meanAUCs[3])
DAUC
