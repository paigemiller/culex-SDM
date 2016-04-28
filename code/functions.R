#Script with support functions of Anopheles arabiensis modeling project
#Created: May 5, 2011 by John M. Drake (jdrake@uga.edu)

#function to get species coordinates
get.sp.coords <- function(sp, plot=FALSE, bbox=NULL){
  require(sp)
  data0<-read.csv('AnGamb_spp.csv')
  data<-data.frame(longitude=data0$LONG[data0$SPECIES==sp], latitude=data0$LAT[data0$SPECIES==sp])
  if(plot==TRUE) plot(data)
  return(SpatialPoints(data, proj4string=CRS('+proj=longlat'), bbox=bbox))
}

#function to rescale a raster by subtracting mean and dividing by standard deviation
rescale.raster <- function(data0, data1=NULL){
  # Rescale a raster by subtracting mean and dividing by standard deviation
  # 
  # Args:
  #   data0: input raster to learn scaling
  #   data1 (optional): second raster to rescale according to the distribution of data0
  #
  # Returns:
  #   A rescaled raster the size of the input
  if(missing(data1)){
    raster <- (data0-cellStats(data0,'mean'))/cellStats(data0,'sd')
  }else{
    raster <- (data1-cellStats(data0,'mean'))/cellStats(data0,'sd')
  }
  return(raster)  
}

#function to rescale a raster by empirical cumulative distribution function
raster.ecdf <- function(data0, data1=NULL){
  # Compute empirical cumulative distribution function on values in a raster
  # Applies ecdf to raster (optionally a second raster) to return a scaled version
  #
  # Args:
  #   data0: input raster to learn scaling
  #   data1 (optional): second raster to rescale according to the ecdf of data0
  #
  # Returns:
  #   A rescaled raster the size of the input
  vals0 <- getValues(data0)
  fun <- ecdf(vals0)
  if(missing(data1)){
    raster <- raster(data0)
    vals <- vals0
  }else{
    raster <- raster(data1)
    vals <- getValues(data1)
  }
  values(raster) <- fun(vals)
  return(raster)
}

#function to take a list of coordinates and buffer radius and return a polygon of buffered area, plot optiona
buffer.to.polygon <- function(coordinates, radius, plot=FALSE){
  require(spatstat)
  c <- coordinates  #local variable for coordinates
  rad <- radius     #local variable for radius
  #now create disc polygons
  polys<-list() 
  for(i in 1:nrow(c)) {
    discbuff<-disc(radius=rad, centre=c(c$longitude[i], c$latitude[i])) 
    discpoly<-Polygon(rbind(cbind(discbuff$bdry[[1]]$x, y=discbuff$bdry[[1]]$y), c(discbuff$bdry[[1]]$x[1], y=discbuff$bdry[[1]]$y[1])))
    polys<-c(polys, discpoly)
  }
  #now list discs
  spolys<-list() 
  for(i in 1:length(polys)) { 
    spolybuff<-Polygons(list(polys[[i]]), ID=row.names(c)[i]) 
    spolys<-c(spolys, spolybuff) 
  } 
  #conver <- raster(data.not.arabiensis)            # create a RasterLayer with the extent of data.arabiensis
  polybuff<-SpatialPolygons(spolys) 
  #plot
  if(plot==TRUE){
    plot(polybuff) 
    points(c, col="red") 
  }
  return(polybuff)
}

#function to thin points from an underdispersed (cluterest) dataset to a given resolution, plot optional
thin.points <- function(points, resolution, plot=TRUE) {
  r <- raster(points)             # create a RasterLayer with the extent of points
  res(r) <- resolution            # set the resolution of the cells
  cell <- cellFromXY(r, points)   # get the cell number for each point
  dup <- duplicated(cell)         # find duplicates
  data.thin <- points[!dup, ]     # select points from unique cells
  if(plot==TRUE){
    p <- rasterToPolygons(r)      # display the results
    plot(p, border='gray')
    points(points)
    points(data.thin, cex=0.5, col='red', pch='x')
  }
  return(data.thin)
}

lobag <- function(p, a, n.votes=100, C=100, sigma=NULL){
  # Applies the lobag algorithm to data on environmental covariates at presence
  # and absence points
  # Requires 'kernlab' library
  
  # Args:
  #   p: covariate data at presence/occurrence points
  #   a: covariate data at absence/background points
  #   n.votes: Number of bootstraps
  #   C: tuning parameter of C-svm
  #   sigma: optional tuning parameter of rbf kernel; uses default if not supplied
  
  # Returns:
  #   A list of length n.votes containing the trained models 
  
  require(kernlab)
  labels <- c(rep(1, dim(p)[1]), rep(0, dim(a)[1]))
  models <- list()
  for(i in 1:n.votes){
    data <- rbind(p[sample.int(n=dim(p)[1], replace=TRUE),], a[sample.int(n=dim(a)[1],
                                                                          replace=TRUE),])
    if(!is.null(sigma)){
      models[[i]] <- ksvm(labels~., data = cbind(data, labels), type='C-svc', C=C,
                          kpar=list(sigma=sigma))
    }
    else{
      models[[i]] <- ksvm(labels~., data = cbind(data, labels), type='C-svc', C=C)
    }
  }
  return(models)
}

lobag.oc <- function(p, n.votes=100, nu=0.01, sigma=NULL){
  # Applies the lobag algorithm to data on environmental covariates at presence points
  # Requires 'kernlab' library
  
  # Args:
  #   p: covariate data at presence/occurrence points
  #   n.votes: Number of bootstraps
  #   nu: tuning parameter of nu-svm
  #   sigma: optional tuning parameter of rbf kernel; uses default if not supplied
  
  # Returns:
  #   A list of length n.votes containing the trained models
  
  require(kernlab)
  labels <- c(rep(1, dim(p)[1]))
  models <- list()
  for(i in 1:n.votes){
    data <- rbind(p[sample.int(n=dim(p)[1], replace=TRUE),])
    if(!is.null(sigma)){
      models[[i]] <- ksvm(labels~., data = cbind(data, labels), type='one-svc', nu=nu,
                          kpar=list(sigma=sigma))
    }
    else{
      models[[i]] <- ksvm(labels~., data = cbind(data, labels), type='one-svc', nu=nu)
    }
  }
  return(models)
}

tune.lobag.oc <- function(p, a, n.votes=512){
  # Uses 10-fold cross-validation to tune the lobag-oc algorithm to data on environmental
  # covariates at presence points
  # Requires 'kernlab' library
  
  # Args:
  #   p: covariate data at presence/occurrence points
  #   a: covariate data at documented absence points
  #   n.votes: Number of bootstraps
  
  # Returns:
  #   A list of length n.votes containing the trained models 
  
  require(kernlab)
  k <- 10
  max.auc <- 0
  n <- dim(p)[1]
  bin <- floor(n/k)
  for(sigma in 2^seq(-2, -9, by=-1)){
    for(nu in 2^seq(-2,-9, by=-1)){
      output.lobag <- list()      
      for(i in 0:(k-1)){
        ids <- (i*bin+1):((i+1)*bin)
        train <- p[-ids,]
        test1 <- p[ids,]
        test0 <- a        
        labels <- c(rep(1, dim(train)[1]))
        models <- list()        
        # fit ensemble
        for(j in 1:n.votes){
          data <- train[sample.int(n=dim(train)[1], replace=TRUE),]
          models[[j]] <- ksvm(labels~., data = cbind(data, labels), type='one-svc', nu=nu,
                              kpar=list(sigma=sigma))
        }       
        # generate predictions and evaluate
        test.lobag <- predictSvm(p=test1, a=test0, model=models) 
        output.lobag[[(i+1)]] <- evaluate(p=as.numeric(test.lobag$p.out),
                                          a=as.numeric(test.lobag$a.out))
      }     
      mean.auc<-mean(unlist(lapply(output.lobag, function(x) x@auc)))      
      if(mean.auc > max.auc){
        max.auc <- mean.auc
        best.model <- models
        best.nu <- nu
        best.sigma <- sigma
      }
    }    
  }
  return(list(auc=max.auc, model=best.model, nu=best.nu, sigma=best.sigma))
}

predictSvm <- function(p, a=NULL, model){
  # Applies a set of base models from lobag function to new data and returns the
  # model prediction
  # Requires 'kernlab' library
  
  # Args:
  #   p: covariate data at primary (presence) points to evaluate
  #   a: covariate data at absence points (optional)
  #   model: output of lobag
  
  # Returns:
  #   p.out: a list of numerical scores at presence points
  #   a.out: a list of numerical scores at absence points
  
  require(kernlab)
  n.votes <- length(model)
  predicts.presence <- p[,1]*0 
  predicts.absence <- NULL
  if(!is.null(a)) predicts.absence <- a[,1]*0
  for(i in 1:n.votes){
    if(!is.null(predicts.presence)) predicts.presence <- predicts.presence + predict(model[[i]], p)
    if(!is.null(predicts.absence)) predicts.absence <- predicts.absence + predict(model[[i]], a)
  }
  return(list(p.out = predicts.presence/n.votes, a.out = predicts.absence/n.votes))
}


##Distance bassed environmental
filterByProximityObsolete <- function(xy, dist, mapUnits = F) {
  #xy can be either a SpatialPoints or SPDF object, or a matrix
  #dist is in km if mapUnits=F, in mapUnits otherwise
  if (!mapUnits) {
    d <- spDists(xy,longlat=F)
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

B.A.Lobag <- function(shittyDataStructure, n.votes){ #dumb function that makes things only slightly better looking
  #makes the data structure slightly less shitty and runs LoBag on species data at a particular distance
  shitList <- list()
  for (j in 1:nspp){
    shitList[[j]] <- lobag.oc(shittyDataStructure[[j]], n.votes=250)
  }
  return(shitList)
}

# getAnswers <- function(nspp, d, pca){ #function to gen spp, thin data, run models, calc & return AUCs
#   #nspp=2   #just for now ##################
#   #nPts <- 2000 #number of sample points
#   #d <- distances
#   
#   ##generates nspp
#   spList <- as.list(seq(1,nspp))
#   sp <- lapply(X=spList, FUN=generateRandomSp, raster.stack=environmental.data.rs, approach="pca", plot=FALSE) #doesn't account for random seeds... is that bad
#   
#   ##sample presence & absence points
#   set.seed(124)
#   sampTrain <- lapply(sp, FUN=sampleOccurrences, n=nPts) #training ;points
#   sampTest <- lapply(sp, FUN=sampleOccurrences, n=nPts, type="presence-absence") #testing points
#   sampTrain <- lapply(sampTrain, FUN=function(x){x[[1]][,1:2]}) #only xy coordinates for training data
#   sampTest <- lapply(sampTest, FUN=function(x){x[[1]][,1:3]}) #xy and presence points for testing data
#   
#   ##extract the environmental variables for the training and testing points
#   train.x <- lapply(sampTrain, FUN=extract, x=environmental.data.rs) #environmental data at training points
#   test.x <- lapply(sampTest, FUN=function(y){extract(y[,1:2], x=environmental.data.rs)})
#   
#   ##Establish 3 points (roughly equivalent to the three clusters in our mosquito data) and calculate distances of each training point to each of them
#   Pt1<-c(30,30)
#   Pt2<-c(10,8)
#   Pt3<-c(35,5)
#   Dist1 <- lapply(sampTrain, FUN=function(x){apply(x, 1, function(x){dist(rbind(x, Pt1))})})
#   Dist2 <- lapply(sampTrain, FUN=function(x){apply(x, 1, function(x){dist(rbind(x, Pt2))})})
#   Dist3 <- lapply(sampTrain, FUN=function(x){apply(x, 1, function(x){dist(rbind(x, Pt3))})})
#   
#   ## Pick the smallest distance to a chosen point for each of the training points
#   speDist <- list()
#   for (i in 1:nspp){
#     speDist[[i]] <- cbind(Dist1[[i]], Dist2[[i]], Dist3[[i]])
#   }
#   DistMin <- lapply(speDist, FUN=function(x){apply(x, 1, min)})
#   
#   ##sample from the training points randomly and according to the minimum distance to a chosen point, get environmental variables at those points
#   set.seed(345)
#   Rand <- lapply(DistMin, FUN=function(x){sample(1:length(x), size=round(nPts/4), replace=F)})
#   Samp <- lapply(DistMin, FUN=function(x){sample(1:length(x), size=round(nPts/4), replace=F, prob=1/((x)^30))})
#   RandPts <- list(); SampPts <- list() 
#   for (i in 1:nspp){
#     RandPts[[i]] <- sampTrain[[i]][Rand[[i]], ]
#     SampPts[[i]] <- sampTrain[[i]][Samp[[i]], ]
#   }
#   RandX <- list() ; SampX <- list()
#   #environmental data
#   for (i in 1:nspp){
#     RandX[[i]] <- train.x[[i]][Rand[[i]], ] #environmental data at random points
#     SampX[[i]] <- test.x[[i]][Samp[[i]], ] #environmental data at biased sample points
#   }
#   
#   ##Thinned points
#   uniquePCA <- lapply(SampX, FUN=function(x){cbind(ID=1:length(x[[1]]), PCA=stats::predict(pca, x))})
#   ThinX <- list() ; envD <- list()
#   for (i in 1:length(d)){
#     thin <- lapply(uniquePCA, FUN=function(x){filterByProximity(x[,2:3], dist=d[i])}) #list of length(d) where each list item is of lenght nspp
#     for (j in 1:nspp){
#       env <- uniquePCA[[j]][,2:3]
#       envD[[j]] <- SampX[[j]][which(env[,1] %in% thin[[j]][[1]][,1]),] #environmental data for each species j at distance D
#     }
#     ThinX[[i]] <- envD
#   }
#   
#   ##train models on random, biased samples, biased corrections of samples, and thinned samples
#   LOB.FULL <- lapply(train.x, FUN=lobag.oc, n.votes=250) 
#   LOB.RAND <- lapply(RandX, FUN=lobag.oc, n.votes=250) 
#   LOB.BIAS <- lapply(SampX, FUN=lobag.oc, n.votes=250)
#   LOB.THIN <- list() ; dat <- list()
#   for (i in 1:length(d)){
#     dist <- ThinX[[i]]
#     for (j in 1:nspp){
#       dat[[j]]<-lobag.oc(dist[[j]], n.votes=250)
#     }
#     LOB.THIN[[i]] <- dat
#   }
#   #LOB.THIN <- lapply(ThinX, FUN=function(y){lapply(y, FUN=lobag.oc, n.votes=250)}) #couldn't make this work :(
#   
#   ##get vals for all the models
#   AUC.LOB.FULL<- list(); AUC.LOB.RAND <- list(); AUC.LOB.BIAS <- list()
#   for (i in 1:nspp){
#     LOB <- predictSvm(model=LOB.FULL[[i]],test.x[[i]]) #list of predictions
#     pred <- ROCR::prediction(LOB$p.out, sampTest[[i]][,3]) #list of prediction items 
#     AUC.LOB.FULL[[i]] <- unlist(ROCR::performance(pred, "auc")@y.values) 
#     
#     LOB <- predictSvm(model=LOB.RAND[[i]], test.x[[i]])
#     pred <- ROCR::prediction(LOB$p.out, sampTest[[i]][,3]) #list of prediction items 
#     AUC.LOB.RAND[[i]] <- unlist(ROCR::performance(pred, "auc")@y.values) #list of AUC values
#     
#     LOB <- predictSvm(model=LOB.BIAS[[i]], test.x[[i]])
#     pred <- ROCR::prediction(LOB$p.out, sampTest[[i]][,3]) #list of prediction items 
#     AUC.LOB.BIAS[[i]] <- unlist(ROCR::performance(pred, "auc")@y.values) #list of AUC values
#   }
#   #AUC vals for env dist thinning
#   dat <- list() ; LOB <- list() ; AUC.LOB.THIN <- list()
#   for (i in 1:length(d)){
#     modD <- LOB.THIN[[i]]
#     for (j in 1:nspp){
#       spp <- modD[[j]]
#       LOB[[j]] <- predictSvm(model=spp, test.x[[j]])
#       pred <- ROCR::prediction(LOB[[j]]$p.out, sampTest[[j]][,3]) #list of prediction items 
#       dat[[j]] <- unlist(ROCR::performance(pred, "auc")@y.values) #list of AUC values
#     }
#     AUC.LOB.THIN[[i]] <- dat
#   }
#   return(list(AUC.FULL, AUC.RAND, AUC.BIAS, AUC.THIN))
# }