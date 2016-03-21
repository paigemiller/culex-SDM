library(geosphere)
library(spThin)
library(dismo)
library(ROCR)
library(foreach)
library(doParallel)

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
load('C:/Users/Robert/OneDrive/Documents/Robbie/8910ENM/Project/swiss-vegetation.Rdata')
x[,2]<-x[,2]/10000
x[,3]<-x[,3]/10000

trainsample<-sample(c(1:550),450)
xtrain<-x[x$ID %in% trainsample,]
ytrain<-y[y$ID %in% trainsample,]
xtest<-x[!(x$ID %in% trainsample),]
ytest<-y[!(x$ID %in% trainsample),]

AUCOrig<-rep(NA, 106)
AUCBiased<-rep(NA, 106)
AUCThinned<-rep(NA, 106)
AUCResamp<-rep(NA, 106)
AUCRand<-rep(NA,106)
AUCOrigMax<-rep(NA, 106)
AUCBiasedMax<-rep(NA, 106)
AUCThinnedMax<-rep(NA, 106)
AUCResampMax<-rep(NA, 106)
AUCRandMax<-rep(NA,106)
for (z in c(42,8,70,3,51,14,63,68)){
#foreach (z=2:107) %do% {
  xtrainPl<-xtrain[ytrain[,z]==1,]

  set.seed(333)
  Point1<-c(runif(1, 55.5,58.5),runif(1, 11.5,15.5))
  Dist1<-rep(NA, nrow(xtrainPl))
  for (i in 1:nrow(xtrainPl))
  {
    Dist1[i]<-euc.dist(xtrainPl[i,2:3], Point1)
  }
  set.seed(999)
  Point2<-c(runif(1, 55.5,58.5),runif(1, 11.5,15.5))
  Dist2<-rep(NA,nrow(xtrainPl))
  for (i in 1:nrow(xtrainPl))
  {
    Dist2[i]<-euc.dist(xtrainPl[i,2:3], Point2)
  }


  DistM<-cbind(Dist1,Dist2)
  Distmin<-data.frame(apply(DistM,1,min))
  Distmin<-cbind(xtrainPl$ID, Distmin)


  colnames(Distmin)<-c('ID','Distance')
  set.seed(345)
  
  Rand<-sample(1:nrow(Distmin), size=round(nrow(Distmin))/2, replace=F)
  Sample<-sample(1:nrow(Distmin),size=round(nrow(Distmin)/2), replace=F, prob=1/Distmin$Distance)

  RandX<-xtrainPl[Rand,]
  SampleX<-xtrainPl[Sample,]


  #plot(SampleX[,2], SampleX[,3], xlim=c(55,59), ylim=c(11,16))
  #points(Point1[1], Point1[2],col='red')
  #points(Point2[1], Point2[2],col='red')

  SPEC<-rep(1,nrow(SampleX))

  thinned<-thin(cbind(SampleX,SPEC), lat.col="X", long.col="Y", spec.col='SPEC', thin.par=10, reps=100, locs.thinned.list.return = T, write.files = F, verbose=F)
  set.seed(123)
  resampled<-sample(1:nrow(SampleX), size=nrow(SampleX), replace=T, prob=(Distmin[Distmin$ID %in% SampleX$ID,2]))
  resampled<-SampleX[resampled,]
  ThinnedSample<-SampleX[as.numeric(rownames(thinned[[50]])),]

  background<-xtrain[!(xtrain$ID %in% xtrainPl$ID),]
  MAX<-maxent(x=rbind(xtrainPl[,4:13],background[,4:13]), p=c(rep(1,nrow(xtrainPl)),rep(0,nrow(background))))
  MaxPred<-predict(MAX,xtest[,4:13])
  Pred<-prediction(MaxPred,ytest[,z])
  AUCOrigMax[z-1]<-as.numeric(performance(Pred,'auc')@y.values)
  
  LOB<-lobag.oc(xtrainPl[,4:13],n.votes=250)
  Predictions<-predictSvm(LOB, xtest[,4:13])
  Pred<-prediction(Predictions$p.out, ytest[,z])
  AUCOrig[z-1]<-as.numeric(performance(Pred,'auc')@y.values)
  
  background<-xtrain[!(xtrain$ID %in% RandX$ID),]
  MAX<-maxent(x=rbind(RandX[,4:13],background[,4:13]), p=c(rep(1,nrow(RandX)),rep(0,nrow(background))))
  MaxPred<-predict(MAX,xtest[,4:13])
  Pred<-prediction(MaxPred,ytest[,z])
  AUCRandMax[z-1]<-as.numeric(performance(Pred,'auc')@y.values)
  
  
  LOB<-lobag.oc(RandX[,4:13],n.votes=250)
  Predictions<-predictSvm(LOB, xtest[,4:13])
  Pred<-prediction(Predictions$p.out, ytest[,z])
  AUCRand[z-1]<-as.numeric(performance(Pred,'auc')@y.values)

  background<-xtrain[!(xtrain$ID %in% SampleX$ID),]
  MAX<-maxent(x=rbind(SampleX[,4:13],background[,4:13]), p=c(rep(1,nrow(SampleX)),rep(0,nrow(background))))
  MaxPred<-predict(MAX,xtest[,4:13])
  Pred<-prediction(MaxPred,ytest[,z])
  AUCBiasedMax[z-1]<-as.numeric(performance(Pred,'auc')@y.values)
  
  LOB<-lobag.oc(SampleX[,4:13],n.votes=250)
  Predictions<-predictSvm(LOB, xtest[,4:13])
  Pred<-prediction(Predictions$p.out, ytest[,z])
  AUCBiased[z-1]<-as.numeric(performance(Pred,'auc')@y.values)

  background<-xtrain[!(xtrain$ID %in% ThinnedSample$ID),]
  MAX<-maxent(x=rbind(ThinnedSample[,4:13],background[,4:13]), p=c(rep(1,nrow(ThinnedSample)),rep(0,nrow(background))))
  MaxPred<-predict(MAX,xtest[,4:13])
  Pred<-prediction(MaxPred,ytest[,z])
  AUCThinnedMax[z-1]<-as.numeric(performance(Pred,'auc')@y.values)
  
  LOB<-lobag.oc(ThinnedSample[,4:13],n.votes=250)
  Predictions<-predictSvm(LOB, xtest[,4:13])
  Predictions
  Pred<-prediction(Predictions$p.out, ytest[,z])
  AUCThinned[z-1]<-as.numeric(performance(Pred,'auc')@y.values)

  background<-xtrain[!(xtrain$ID %in% resampled$ID),]
  MAX<-maxent(x=rbind(resampled[,4:13],background[,4:13]), p=c(rep(1,nrow(resampled)),rep(0,nrow(background))))
  MaxPred<-predict(MAX,xtest[,4:13])
  Pred<-prediction(MaxPred,ytest[,z])
  AUCResampMax[z-1]<-as.numeric(performance(Pred,'auc')@y.values)
  
  LOB<-lobag.oc(resampled[,4:13],n.votes=250)
  Predictions<-predictSvm(LOB, xtest[,4:13])
  Pred<-prediction(Predictions$p.out, ytest[,z])
  AUCResamp[z-1]<-as.numeric(performance(Pred,'auc')@y.values)
  
#ThinnedSample
#plot(ThinnedSample[,2], ThinnedSample[,3],  xlim=c(55,59), ylim=c(11,16))


#plot(xtrainPl3[,2], xtrainPl3[,3], xlim=c(55,59), ylim=c(11,16))
#points(Point1[1], Point1[2],col='red')
#points(Point2[1], Point2[2],col='red')
}

AUCMat<-cbind(AUCOrig,AUCBiased,AUCThinned, AUCResamp,AUCRand)
AUCMat.narm<-AUCMat[which(!is.na(AUCMat[,4])),]
mean(AUCMat.narm[,1]-AUCMat.narm[,2])
mean(AUCMat.narm[,1]-AUCMat.narm[,3])
mean(AUCMat.narm[,2]-AUCMat.narm[,3])



AUCMat.final<-AUCMat.narm[which(AUCMat.narm[,1]>AUCMat.narm[,2]),]


DAUCThin<-(AUCMat.final[,3]-AUCMat.final[,2])/(AUCMat.final[,1]-AUCMat.final[,2])
DAUCResamp<-(AUCMat.final[,4]-AUCMat.final[,2])/(AUCMat.final[,1]-AUCMat.final[,2])


mean(AUCMat.final[,1]-AUCMat.final[,2])
mean(AUCMat.final[,1]-AUCMat.final[,3])
mean(AUCMat.final[,2]-AUCMat.final[,3])
mean(AUCMat.final[,2]-AUCMat.final[,4])

AUCMatMax<-cbind(AUCOrigMax,AUCBiasedMax,AUCThinnedMax, AUCResampMax,AUCRandMax)
AUCMatMax.narm<-AUCMatMax[which(!is.na(AUCMatMax[,4])),]
AUCMatMax.final<-AUCMatMax.narm[which(AUCMatMax.narm[,1]>AUCMatMax.narm[,2]),]

mean(AUCMatMax.final[,1]-AUCMatMax.final[,2])
mean(AUCMatMax.final[,1]-AUCMatMax.final[,3])
mean(AUCMatMax.final[,2]-AUCMatMax.final[,3])
mean(AUCMatMax.final[,2]-AUCMatMax.final[,4])
