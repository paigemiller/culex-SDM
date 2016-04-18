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
filterByProximity <- function(xy, dist, mapUnits = F) {
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

