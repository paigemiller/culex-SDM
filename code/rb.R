rb <- function(x, v=100, d=1, p=0.5){
  # Basic version of the range bagging algorithm 
  # Requires geometry package
  #
  # Args:
  #   x: covariate data at observations
  #   v: number of votes
  #   d: dimension of ranges to bag
  #
  # Returns:
  #   A list of length v containing the individual models 
  require(geometry)
  models <- list()
  n <- dim(x)
  for(i in 1:v){
    vars <- sample.int(n[2], size=d, replace=FALSE)
    x0 <- x[,vars]
    
    if(d==1){
      x1 <- x0[sample(n[1],ceiling(p*n[1]), replace=FALSE)]
      models[[i]] <- list(vars=vars, endpoints=c(min(x1), max(x1)), data=x1)
    }
    else{
      x1 <- x0[sample(n[1],ceiling(p*n[1]), replace=FALSE),]
      #THIS DOESNT REALLY DO ANYTHING USEFUL SINCE WE REFIT CONVULL IN FXN RB
      idx <- unique(as.vector(convhulln(x1, options='Pp')))
      endpoints <- x1[idx,]
      models[[i]] <- list(vars=vars, endpoints=endpoints, data=unique(x1))
    }    
  }
  return(models)
}


rb.test <- function(models, x.new){
  # Test function for a point to determine the fraction of bags in which the point falls
  # Requires geometry package
  #
  # Args:
  #   models: a list of models returned by function rb
  #   x.new: covariate data at new points
  #
  # Returns:
  #   A vector of the same length as x.new containing the fraction of models for each point in which it is included
  v <- length(models)
  d <- ifelse(is.null(dim(models[[1]]$endpoints)), 1, dim(models[[1]]$endpoints)[2])
  n <- dim(x.new)
  out <- numeric(n[1])
  for(i in 1:v){
    #    print(i) # counter for troubleshooting
    if(d==1){
      test.pts <- (models[[i]]$endpoints[1] < x.new[,models[[i]]$vars]) & (x.new[,models[[i]]$vars] < models[[i]]$endpoints[2])
      out <- out + test.pts      
    }else{
      test.dat <- as.matrix(x.new[,models[[i]]$vars])
      tri.pts <- tsearchn(as.matrix(models[[i]]$data), delaunayn(models[[i]]$data), test.dat)
      #tri.pts <- tsearchn(as.matrix(models[[i]]$endpoints), delaunayn(models[[i]]$endpoints), test.dat)
      test.pts <- !is.na(tri.pts$p[,1])
      out <- out + test.pts      
    }
  }
  return(out/v)
}