# generate raster index
rast_index<-function(indx,geom=c("lon","lat"),resolution=1,crs="+proj=longlat +datum=WGS84",extent=1,field=""){
  library(terra)
  pts <- vect(indx, geom=geom, crs=crs)
  r<-rast()
  ext(r)<-ext(pts)
  res(r)<-resolution
  r<-extend(r,ext(r)+(extent*resolution))
  pir<-rasterize(pts,r,field=field)
}
# get the mean of duplicated points
dup_mean<-function(ind.rast, indx){
  library(terra)
  pt<-indx[,c("lon","lat")]
  ll<-cbind(cellFromXY(ind.rast,pt),indx)
  dp<-unique(ll[duplicated(ll[,1]),1])
  if(length(dp)>0){
    wodp<-ll[!duplicated(ll[,1]),]
    for(i in seq_along(dp)){
      wodp[wodp[,1]==dp[i],4]<-mean(ll[ll[,1]==dp[i],4])
    }
  } else {
    wodp<-indx
  }
  return(wodp[,2:4])
}