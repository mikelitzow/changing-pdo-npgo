#plot_geo_map.r
#Plots geographic maps of codes for each node, called in parent script runSOM2_xyf.r

#library('RColorBrewer')
#source('filled_contour3.R')
#source('filled_legend.R')

library(fields)  #for image.plot to include legend, for two.colors (or tim.colors)
library(maps)
library(mapdata)
library(pracma)   #for meshgrid function

# Get codes for each variable
if (compsize == 1) {
  codes_out1 <- codes_out
} 
if (compsize == 2 ) {
  codes_out1 <- codes_out[[1]]
  codes_out2 <- codes_out[[2]]
} 
if (compsize == 3) {
  codes_out1 <- codes_out[[1]]
  codes_out2 <- codes_out[[2]]
  codes_out3 <- codes_out[[3]]
}
if (compsize == 4) {
  codes_out1 <- codes_out[[1]]
  codes_out2 <- codes_out[[2]]
  codes_out3 <- codes_out[[3]]
  codes_out4 <- codes_out[[4]]
}

# Make plots for each variable (if the variable isn't spatially explicit, this won't make sense, but that's okay) - try later with ggplot??
zmax <- max(abs(min(codes_out1)),max(codes_out1))
zmin <- -zmax
dev.new()
par(mfrow=c(sdim2,sdim1),
    oma = c(5,1,0,0) + 0.1,
    mar = c(0,2,1,1) + 0.1)
MG1 <- meshgrid(lonvecM1, latvecM1)
for (node in c(1:nodes)){
  dummyT <- matrix(data=NA,nrow=length(latvecM1),ncol=length(lonvecM1))
  dummyT[seaindsM1]<-codes_out1[node,]
  dummy <- t(dummyT)
  image.plot(x=MG1$X, y=MG1$Y, z=dummyT, xlim=c(lonmin,lonmax), ylim=c(latmin,latmax), zlim=c(zmin,zmax), col=two.colors(start='blue', end = 'red', middle='white'), main=paste(sprintf('%s', grids[1]), 'node', node))
  contour(x=lonvecM1, y=latvecM1, z=dummy, zlim=c(zmin,zmax), add=T)
  map('world2Hires', fill=F, add=T, lwd=2)
  }
par(new = "TRUE",plt = c(0.85,0.9,0.25,0.85),las = 1,cex.axis = 1)

if (compsize > 1) {
  zmax <- max(abs(min(codes_out2)),max(codes_out2))
  zmin <- -zmax
  dev.new()
  par(mfrow=c(sdim2,sdim1),
      oma = c(5,1,0,0) + 0.1,
      mar = c(0,2,1,1) + 0.1)
MG2 <- meshgrid(lonvecM2, latvecM2)
  for (node in c(1:nodes)){
    dummyT <- matrix(data=NA,nrow=length(latvecM2),ncol=length(lonvecM2))
    dummyT[seaindsM2]<-codes_out2[node,]
    dummy <- t(dummyT)
    image.plot(x=MG2$X, y=MG2$Y, z=dummyT, xlim=c(lonmin,lonmax), ylim=c(latmin,latmax), zlim=c(zmin,zmax), col=two.colors(start='blue', end = 'red', middle='white'), main=paste(sprintf('%s', grids[2]), 'node', node))
    contour(x=lonvecM2, y=latvecM2, z=dummy, zlim=c(zmin,zmax), add=T)
    map('world2Hires', fill=F, add=T, lwd=2)
  }
  par(new = "TRUE",plt = c(0.85,0.9,0.25,0.85),las = 1,cex.axis = 1)
}

if (compsize > 2) {
  zmax <- max(abs(min(codes_out3)),max(codes_out3))
  zmin <- -zmax
  dev.new()
  par(mfrow=c(sdim2,sdim1),
      oma = c(5,1,0,0) + 0.1,
      mar = c(0,2,1,1) + 0.1)
MG3 <- meshgrid(lonvecM3, latvecM3)
  for (node in c(1:nodes)){
    dummyT <- matrix(data=NA,nrow=length(latvecM3),ncol=length(lonvecM3))
    dummyT[seaindsM3]<-codes_out3[node,]
    dummy <- t(dummyT)
    image.plot(x=MG3$X, y=MG3$Y, z=dummyT, xlim=c(lonmin,lonmax), ylim=c(latmin,latmax), zlim=c(zmin,zmax), col=two.colors(start='blue', end = 'red', middle='white'), main=paste(sprintf('%s', grids[3]), 'node', node))
    contour(x=lonvecM3, y=latvecM3, z=dummy, zlim=c(zmin,zmax), add=T)
    map('world2Hires', fill=F, add=T, lwd=2)
  }
  par(new = "TRUE",plt = c(0.85,0.9,0.25,0.85),las = 1,cex.axis = 1)
}

if (compsize > 3) {
  zmax <- max(abs(min(codes_out4)),max(codes_out4))
  zmin <- -zmax
  dev.new()
  par(mfrow=c(sdim2,sdim1),
      oma = c(5,1,0,0) + 0.1,
      mar = c(0,2,1,1) + 0.1)
MG4 <- meshgrid(lonvecM4, latvecM4)
  for (node in c(1:nodes)){
    dummyT <- matrix(data=NA,nrow=length(latvecM4),ncol=length(lonvecM4))
    dummyT[seaindsM4]<-codes_out4[node,]
    dummy <- t(dummyT)
    image.plot(x=MG4$X, y=MG4$Y, z=dummyT, xlim=c(lonmin,lonmax), ylim=c(latmin,latmax), zlim=c(zmin,zmax), col=two.colors(start='blue', end = 'red', middle='white'), main=paste(sprintf('%s', grids[4]), 'node', node))
    contour(x=lonvecM4, y=latvecM4, z=dummy, zlim=c(zmin,zmax), add=T)
    map('world2Hires', fill=F, add=T, lwd=2)
    }
  par(new = "TRUE",plt = c(0.85,0.9,0.25,0.85),las = 1,cex.axis = 1)
}
