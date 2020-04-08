## plot_geo_map.r
## Emily Norton, 8/7/18
## Plots geographic maps of codes for each node, called in parent script runSOM_xyf_super.r

library(fields)  #for image.plot to include legend, for two.colors (or tim.colors)
library(maps)
library(mapdata)
library(pracma)   #for meshgrid function

for (d in 1) {
	
	# Get Codes
	codestemp <- codes_out[[d]]
	
	# Make plots for each variable (if the variable isn't spatially explicit, this won't make sense, but that's okay) 
	zmax <- max(abs(min(codestemp)),max(codestemp))
	zmin <- -zmax
	dev.new()
	par(mfrow=c(sdim2,sdim1),
    	oma = c(5,1,0,0) + 0.1,
    	mar = c(0,2,1,1) + 0.1)
	MG <- meshgrid(lslonvec[[d]], lslatvec[[d]])
	
	for (node in c(1:nodes)){
  		dummyT <- matrix(data=NA,nrow=length(lslatvec[[d]]),ncol=length(lslonvec[[d]]))
  		dummyT[lsseainds[[d]]]<-codestemp[node,]
  		dummy <- t(dummyT)
  	
  		image.plot(x=MG$X, y=MG$Y, z=dummyT, xlim=c(lonmin,lonmax), ylim=c(latmin,latmax),
  			col=two.colors(start='blue', end = 'red', middle='white'), main=paste(sprintf('%s', grids[[d]]), 'node', node))
  			contour(x=lslonvec[[d]], y=lslatvec[[d]], z=dummy, zlim=c(zmin,zmax), add=T)
  		map('world2Hires', fill=F, add=T, lwd=2)
  	}
	
	par(new = "TRUE",plt = c(0.85,0.9,0.25,0.85),las = 1,cex.axis = 1)

}
