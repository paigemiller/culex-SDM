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
