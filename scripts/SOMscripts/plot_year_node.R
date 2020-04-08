## plot_year_node.r 
## Emily Norton, 8/7/18
## Generates a hovmoller-style plot for plotting the node each year was sorted into

require("grDevices") # for colors

#Load in the data with the maps$unit_classif (i.e. the node number for each year) - generated in runSOM_xyf_super_FINAL.r
maps_uclass <- as.matrix(maps_out$unit.classif) 

# Pre-allocate matrix of m years (rows) and n nodes (columns)
	maps_z <- matrix(data=NA,nrow=length(maps_uclass),ncol=nodes)
	
	for (entry in c(1:length(maps_uclass))){
  	maps_z[entry,maps_uclass[entry]]<-1
	}

# Generate Hovmoller plot
	dev.new()

	image(c(1:nodes), seq(yearmin, yearmax+(1-1/nmon), by=(1/nmon)), t(maps_z),   	# using seq(...) allows the months to be evenly distributed within the year
      ylim = c(yearmax+(1-1/nmon), yearmin),                                  		# add new max
      col = "black",                                                          		# Future: Add different colors for different months
      ylab = "Year", xlab = "Node", main = title_name)
