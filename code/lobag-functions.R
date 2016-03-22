# This file contains code for four functions used to fit and predict species' ecological niches using the lobag and lobag-oc models described in Drake (2014).
# (c) 2012 by John M. Drake
# This file created: 24 April 2014 by John M. Drake (email: jdrake@uga.edu)
# Dependencies: 'kernlab'
# Reference: Drake, J.M. 2014. Ensemble algorithms for ecological niche modeling from presence-background and presence-only data. Ecosphere.

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
