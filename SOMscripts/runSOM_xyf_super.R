## runSOM_xyf_super.r
## Emily Norton, 8/7/18
## This script uses the som, xyf, or supersom function of the kohonen package to create self-organizing maps for one or more variables

#clear any objects
	rm(list = ls(all.names = TRUE))
#clear any potentially interfering loaded packages
	lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

# Set seed (optional)
	set.seed(10)   

# Set number of variables to compare 
	compsize <- 3   #1 = univariable "som", 2 = two co-variable "xyf", 3 or more = multivariable "supersom"

# Set variable filenames to load - only include data for the sea points
# NOTE: these files must contain headers and 'year' column, and must *at least* contain the temporal range of interest
	filein <- NULL
	filein[1] <- 'SOMscripts/InputData/MonthlyAvg11to3_forYears1948to2018_Bounds20to65N_120to255E_skt_notdetrended.csv'
	filein[2] <- 'SOMscripts/InputData/PDO.csv'     # PDO averaged Nov-March
	filein[3] <- 'SOMscripts/InputData/NPGO.csv'     # NPGO averaged Nov-March


# Set variable grids: "HGT_A", "SLP"/"HGT" (both same grid), "SKT", "MLD", "UFLX", "VFLX", "CURL". 
# If non-spatially explicity (e.g. timeseries of biol V1, PDO, or NPGO), a default grid ("SLP") will be used to generate plots, so enter the variable name below 
	vargrid <- NULL
	vargrid[1] <- "SKT"		# grid for var 1
	vargrid[2] <- "PDO" 		# grid for var 2
	vargrid[3] <- "NPGO"		# grid for var 3


# Set temporal range of interest  
	yearmin <- 1951
	yearmax <- 2018

# Set number of months per year to read in (must be the number of rows per year in the input .csv) 
# NOTE: Each file would need this redundancy
	nmon <- 1

# Set spatial boundaries for the map grid  (requires "year" column, and row header)   
	latmin = 20
	latmax = 65
	lonmin = 120
	lonmax = 255

# Set dimensions for the som grid (2D) and shape of somgrid (options: 'rectangular', 'hexagonal')
	sdim1 <- 2 
	sdim2 <- 3
	sshape <- 'rectangular' 

# Set relative weights for each variable 
	wei <- NULL
	wei[1] <- 998   # weight for var 1
	wei[2] <- 1   # weight for var 2
	wei[3] <- 1   # weight for var 3



# Select types of plots to generate (options: 'T' = true or 'F' = false)
	booYearNode <- 'T'		# Hovmoller-style plot of years and nodes
	booGeoMap <- 'T'    	# Spatially-explicit map for each variable and node

# Plot title and name to save plots
	title_name <- "Year Nodes"    # for Year_Node plot
	plot_fname <- "plot_filename"

# Save SOM codes or maps as .csv files? (options: 'T' = true or 'F' = false)
	saveCodes <- 'F'
	saveMaps <- 'F'


##----------Shouldn't need to change much below this line ----------------

#Load kohonen package (load other packages later)
	library(tidyverse)   #for dplyr and ggplot - load first because of conflicts
	library(kohonen)  #for SOMs

#Calculate the number of nodes that all codes are grouped into
	nodes <- sdim1*sdim2 

# Pre-allocate list variables
	grids <- list()
	weights <- NULL
	lsDfS <- list()              
	lsseainds <- list()
	lslatvec <- list()
	lslonvec <- list()


for (d in 1:compsize) {
	
# Set up grid names and weights, for naming files and running som later
	gtemp <- vargrid[d]					
	grids <- c(grids,gtemp)                 
	wtemp <- wei[d]
	weights <- c(weights, wtemp)


#Load csv file(s) and grab only the values of interest for SOM comparison (i.e. not headers)
	datain <- read_csv(filein[d])  
	data_forSOM <- datain %>%   #get rid of year column...may want to prepare in other ways, e.g. make sure all of the same years are used for both files, etc.
  		filter(floor(year) <= yearmax, floor(year) >= yearmin) %>%        # note: use floor(year) since we have decimal years
  		arrange(year)                                                     # put data in ascending order for plotting later

	lsDfS[[d]] <- data_forSOM


# Set grid filenames to load, with only lat/lon for the sea points, and the indexing info for each point - New user will need to update for their specific grid

	if (sum(vargrid[d] == c("HGT_A", "SLP", "SKT", "MLD", "UFLX", "VFLX", "CURL"))>0) {
		seaindsfile <- sprintf('SOMscripts/InputData/%s_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE_incLat.csv', vargrid[d],latmin, latmax, lonmin, lonmax)
		latvecfile <- sprintf('SOMscripts/InputData/%s_lat_vec_Bounds_OrigGrid_20to65N_120to255E_incLat.csv', vargrid[d])
		lonvecfile <- sprintf('SOMscripts/InputData/%s_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv', vargrid[d])
	} else {   
		print('grid not found, using default SLP grid')
		seaindsfile <- sprintf('SOMscripts/InputData/SLP_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE_incLat.csv', latmin, latmax, lonmin, lonmax)
		latvecfile <- 'SOMscripts/InputData/SLP_lat_vec_Bounds_OrigGrid_20to65N_120to255E_incLat.csv'
		lonvecfile <- 'SOMscripts/InputData/SLP_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv'
	}

# Load grid files in as matrices, and the indices ('seainds') indicating the spatial location for each point - NOTE: 'SLP' (sea level pressure)
# and 'HGT' (geopotential height at 200mbar) have the same grid and sea indices
	lsseainds[[d]] <- as.matrix(read.csv(seaindsfile, header = FALSE)) ## NOTE: weird stuff with HGT_A re: sea only vs including land...
	lslatvec[[d]] <-as.matrix(read.csv(latvecfile, header=FALSE))
	lslonvec[[d]] <-as.matrix(read.csv(lonvecfile, header=FALSE))
	
	## A note on input files: Since the Kohonen som function does not allow NAN points to be included in the variable vector, 
	## points outside of the boundaries and/or on land have been removed from the variable input files. To plot these points 
	## in a pretty way, we require the following files:
		# latvec/lonvec files: vector of lat/lon points from the grid that spans at least the entire plotting domain (required 
		#	for the spatially-explicit plotting); latitude must be monotonically increasing 
		# seainds: indices of the points from your variable input file based on the full grid (required for the spatially-explicit plotting)

}

varnames <- paste(c(grids),collapse='_')   # for saving files


#Run SOMs with 1 ("som"), 2 ("xyf"), or 3 or more variables ("supersom")
  	if (compsize == 1) {
  		som_out <- som(scale(lsDfS[[1]][,-1]),grid=somgrid(sdim1,sdim2,sshape))
	} else if (compsize == 2) {
  		som_out <- xyf(X=scale(lsDfS[[1]][,-1]),Y=scale(lsDfS[[2]][,-1]),user.weights=weights,grid=somgrid(sdim1,sdim2,sshape))
	} else {
  		Comb_dat <- list(scale(lsDfS[[1]][,-1]))
  		for (d in 2:compsize) {
  			new_dat <- list(scale(lsDfS[[d]][,-1]))    
  			Comb_dat <- c(Comb_dat, new_dat)
 		}
  		som_out <- supersom(Comb_dat, user.weights=weights,grid=somgrid(sdim1,sdim2,sshape))
	}


# Get codes for the SOM and save, if desired
	codes_out <- getCodes(som_out)
	
	for (d in 1:compsize) { 
		if (saveCodes == 'T') {
			codestemp <- codes_out[[d]]
  			codes_fname <- sprintf("output_codes%i_SOM_%iby%i%sgrid_%ito%iyears_%s.csv", d, sdim1, sdim2, sshape, yearmin, yearmax, vargrid[d])
  			write.table(codestemp, file = codes_fname, row.names=FALSE,col.names=FALSE,sep=',')
  			codestemp <- NULL
  			codes_fname <- NULL
  		}
  	}

comp.dat <- read.csv("/Users/MikeLitzow 1/Documents/R/FATE2 non-som legacy/SOM.map.data.csv", row.names = 1)

plot.dat <- comp.dat[7:8,]

write.csv(plot.dat, "lat.long.for.SOM.plotting.csv", row.names = F)

# identify NAs
replace <- !is.na(comp.dat[1,])
write.csv(replace, "cells.to.replace.csv", row.names = F)

replace <- read.csv("cells.to.replace.csv")
replace <- as.matrix(replace)
replace <- as.vector(replace)


SOM_out <- as.matrix(codes_out[[1]])

lat.long <- read.csv("lat.long.for.SOM.plotting.csv")

nodes <- matrix(NA, 6, ncol(lat.long))
rownames(nodes) <- c("node1", "node2", "node3", "node4", "node5", "node6")
colnames(nodes) <- colnames(lat.long)
plot.SOM <- as.matrix(rbind(nodes, lat.long))

for(i in 1:6){
plot.SOM[i,replace] <- SOM_out[i,] 
}


dat <- as.data.frame(t(plot.SOM))
colnames(dat)[7:8] <- c("lat", "long")

library(oce)
library(maps)
library(mapdata)
library(fields)

# and an alternate layout!
tiff("figs/Fig 3.tiff", 8, 9, units="cm", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 2
l.cex <- 0.6
l.l <- 0.2
tc.l <- -0.2
l.w <- 0.8

new.col <- oceColorsPalette(64)

par(mar=c(0.5,0.5,0.5,0.5),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(3,2), 
    cex.axis=0.8, cex.lab=0.8, oma=c(0, 0, 0.2, 0.1))

node1 <- as.matrix(tapply(dat$node1, list(dat$lat, dat$long), mean))

x <- as.numeric(colnames(node1))
y <- as.numeric(rownames(node1))

z <- dat$node1
z <- t(matrix(z, length(y))) 
image.plot(x, y, z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l,
           legend.width = l.w,
           axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="grey", drawlabels = F, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("a) Node 1", cex=0.5, adj=0)

z <- dat$node2
z <- t(matrix(z, length(y))) 
image.plot(x, y, z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l,
           legend.width = l.w,
           axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="grey", drawlabels = F, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("b) Node 2", cex=0.5, adj=0)

z <- dat$node3
z <- t(matrix(z, length(y))) 
image.plot(x, y, z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l,
           legend.width = l.w,
           axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="grey", drawlabels = F, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("c) Node 3", cex=0.5, adj=0)

z <- dat$node4
z <- t(matrix(z, length(y))) 
image.plot(x, y, z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l,
           legend.width = l.w,
           axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="grey", drawlabels = F, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("d) Node 4", cex=0.5, adj=0)

z <- dat$node5
z <- t(matrix(z, length(y))) 
image.plot(x, y, z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l,
           legend.width = l.w,
           axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="grey", drawlabels = F, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("e) Node 5", cex=0.5, adj=0)

z <- dat$node6
z <- t(matrix(z, length(y))) 
image.plot(x, y, z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l,
           legend.width = l.w,
           axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="grey", drawlabels = F, lwd=0.5)
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
mtext("f) Node 6", cex=0.5, adj=0)

dev.off()




