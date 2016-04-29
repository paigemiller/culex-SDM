############################code for environmental distance thinning########################
#rm(list = ls()) #run before parallel code
source('functions.R') #filterByProximity function & others
load("culex-v3.RData")

##Ripley's L of PCA space
pca.ppp <- as.ppp(pca$x[,1:2], c(-20,20,-20,20)) #rejects 68 points... why? 
pca.lest <- Lest(pca.ppp, correction="Ripley", nlarge=30000, rmax=5)
plot.fv(pca.lest)

##parameters for thinning
N=10 #number of species
distances <- seq(.2, 5, by=.25)

start <- proc.time() #begins elapsed time of the simulation
envThinDat <- getAnswers(10, distances, pca=pca) 
stopt <- proc.time()
elapsed <- stopt - start
print(elapsed) #prints elapsed time to console

getAnswers <- function(index, d, pca){ #function to gen spp, thin data, run models, calc & return AUCs
  #nspp=2   #just for now ##################
  nPts <- 2000 #number of sample points
  #d <- distances
  
  nspp=index
  
  ##generates nspp
  set.seed(index)
  spList <- as.list(seq(1))###,nspp))
  sp <- lapply(X=spList, FUN=generateRandomSp, raster.stack=environmental.data.rs, approach="pca", plot=FALSE) #doesn't account for random seeds... is that bad
  
  ##sample presence & absence points
  set.seed(124)
  sampTrain <- lapply(sp, FUN=sampleOccurrences, n=nPts) #training ;points
  sampTest <- lapply(sp, FUN=sampleOccurrences, n=nPts, type="presence-absence") #testing points
  sampTrain <- lapply(sampTrain, FUN=function(x){x[[1]][,1:2]}) #only xy coordinates for training data
  sampTest <- lapply(sampTest, FUN=function(x){x[[1]][,1:3]}) #xy and presence points for testing data
  
  ##extract the environmental variables for the training and testing points
  train.x <- lapply(sampTrain, FUN=extract, x=environmental.data.rs) #environmental data at training points
  test.x <- lapply(sampTest, FUN=function(y){extract(y[,1:2], x=environmental.data.rs)})
  
  ##Establish 3 points (roughly equivalent to the three clusters in our mosquito data) and calculate distances of each training point to each of them
  Pt1<-c(30,30)
  Pt2<-c(10,8)
  Pt3<-c(35,5)
  Dist1 <- lapply(sampTrain, FUN=function(x){apply(x, 1, function(x){dist(rbind(x, Pt1))})})
  Dist2 <- lapply(sampTrain, FUN=function(x){apply(x, 1, function(x){dist(rbind(x, Pt2))})})
  Dist3 <- lapply(sampTrain, FUN=function(x){apply(x, 1, function(x){dist(rbind(x, Pt3))})})
  
  ## Pick the smallest distance to a chosen point for each of the training points
  speDist <- list()
  for (i in 1:1){
    speDist[[i]] <- cbind(Dist1[[i]], Dist2[[i]], Dist3[[i]])
  }
  DistMin <- lapply(speDist, FUN=function(x){apply(x, 1, min)})
  
  ##sample from the training points randomly and according to the minimum distance to a chosen point, get environmental variables at those points
  set.seed(345)
  Rand <- lapply(DistMin, FUN=function(x){sample(1:length(x), size=round(nPts/4), replace=F)})
  Samp <- lapply(DistMin, FUN=function(x){sample(1:length(x), size=round(nPts/4), replace=F, prob=1/((x+1)^3))})
  RandPts <- list(); SampPts <- list() 
  for (i in 1:1){
    RandPts[[i]] <- sampTrain[[i]][Rand[[i]], ]
    SampPts[[i]] <- sampTrain[[i]][Samp[[i]], ]
  }
  RandX <- list() ; SampX <- list()
  #environmental data
  for (i in 1:1){
    RandX[[i]] <- train.x[[i]][Rand[[i]], ] #environmental data at random points
    SampX[[i]] <- test.x[[i]][Samp[[i]], ] #environmental data at biased sample points
  }
  
  ##Thinned points
  uniquePCA <- lapply(SampX, FUN=function(x){cbind(ID=1:length(x[[1]]), PCA=stats::predict(pca, x))})
  ThinX <- list() ; envD <- list()
  for (i in 1:length(d)){
    thin <- lapply(uniquePCA, FUN=function(x){filterByProximity(x[,2:3], dist=d[i])}) #list of length(d) where each list item is of lenght nspp
    for (j in 1:1){
      env <- uniquePCA[[j]][,2:3]
      envD[[j]] <- SampX[[j]][which(env[,1] %in% thin[[j]][[1]][,1]),] #environmental data for each species j at distance D
    }
    ThinX[[i]] <- envD
  }
  
  ##train models on random, biased samples, biased corrections of samples, and thinned samples
  LOB.FULL <- lapply(train.x, FUN=lobag.oc, n.votes=250) 
  LOB.RAND <- lapply(RandX, FUN=lobag.oc, n.votes=250) 
  LOB.BIAS <- lapply(SampX, FUN=lobag.oc, n.votes=250)
  LOB.THIN <- list() ; dat <- list()
  for (i in 1:length(d)){
    dist <- ThinX[[i]]
    for (j in 1:1){
      dat[[j]]<-lobag.oc(dist[[j]], n.votes=250)
    }
    LOB.THIN[[i]] <- dat
  }
  #LOB.THIN <- lapply(ThinX, FUN=function(y){lapply(y, FUN=lobag.oc, n.votes=250)}) #couldn't make this work :(

  ##get vals for all the models
  AUC.LOB.FULL<- list(); AUC.LOB.RAND <- list(); AUC.LOB.BIAS <- list()
  for (i in 1:1){
    LOB <- predictSvm(model=LOB.FULL[[i]],test.x[[i]]) #list of predictions
    pred <- ROCR::prediction(LOB$p.out, sampTest[[i]][,3]) #list of prediction items 
    AUC.LOB.FULL[[i]] <- unlist(ROCR::performance(pred, "auc")@y.values) 
    
    LOB <- predictSvm(model=LOB.RAND[[i]], test.x[[i]])
    pred <- ROCR::prediction(LOB$p.out, sampTest[[i]][,3]) #list of prediction items 
    AUC.LOB.RAND[[i]] <- unlist(ROCR::performance(pred, "auc")@y.values) #list of AUC values
    
    LOB <- predictSvm(model=LOB.BIAS[[i]], test.x[[i]])
    pred <- ROCR::prediction(LOB$p.out, sampTest[[i]][,3]) #list of prediction items 
    AUC.LOB.BIAS[[i]] <- unlist(ROCR::performance(pred, "auc")@y.values) #list of AUC values
  }
  #AUC vals for env dist thinning
  dat <- list() ; LOB <- matrix(NA, nrow=2000, ncol=length(d)) ; AUC.LOB.THIN <- list()
  for (i in 1:length(d)){
    modD <- LOB.THIN[[i]]
    for(j in 1:1)
    {
      spp <- modD[[j]]
      LOB[,i] <- predictSvm(model=spp, test.x[[j]])$p.out
      pred <- ROCR::prediction(LOB[,i], sampTest[[j]][,3]) #list of prediction items 
      dat[[j]] <- unlist(ROCR::performance(pred, "auc")@y.values) #list of AUC values
      
    }
    AUC.LOB.THIN[[i]] <- dat
  }
  LOBMEAN<-apply(LOB, MARGIN=1, FUN=mean)
  pred<-ROCR::prediction(LOBMEAN, sampTest[[1]][,3])
  AUC.LOB.MEAN<-unlist(ROCR::performance(pred, "auc")@y.values)
  return(list(unlist(AUC.LOB.FULL), unlist(AUC.LOB.RAND), unlist(AUC.LOB.BIAS), unlist(AUC.LOB.THIN), unlist(AUC.LOB.MEAN)))
}

##parallel code

cluster <- makeCluster(10, type="SOCK") # Creates a set of n reps of R running in parallel and communicating over sockets
# load required libraries on all cluster nodes
clusterEvalQ(cluster, library(virtualspecies))
clusterEvalQ(cluster, library(dismo))
clusterEvalQ(cluster, library(ROCR))
clusterEvalQ(cluster, library(spThin))
clusterEvalQ(cluster, library(sp))
clusterEvalQ(cluster, source('functions.R'))
clusterExport(cluster, c("environmental.data.rs"))
clusterExport(cluster, c("pca"))
#clusterExport(cluster, ls())

clusterSetupRNG(cluster) # set up random number generator

##parameters for thinning
N=10 #number of species
distances <- seq(.2, 5, by=.25)

start <- proc.time() #begins elapsed time of the simulation
envThinDat100 <- snow::clusterApplyLB(cluster, fun=getAnswers, 1:100, d=distances, pca=pca) # disperse individuals and spread the disease 
stopt <- proc.time() #ends elapsed time of the simulation
elapsed <- stopt - start
print(elapsed) #prints elapsed time to console
stopCluster(cluster) # terminate cluster

##parameters for thinning
#N=10 #number of species
distances <- seq(.2, 5, by=.25)

start <- proc.time() #begins elapsed time of the simulation
envThinDat <- getAnswers(1, distances, pca=pca) 
stopt <- proc.time()
elapsed <- stopt - start
print(elapsed) #prints elapsed time to console








