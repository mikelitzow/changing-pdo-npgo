#plot_year_node.r 
#Generates a hovmoller-esque plot for plotting the node each year was sorted into

#Load in the data with the maps$unit_classif (i.e. the node number for each year) - generated in mapSOMs.r
maps_uclass <- as.matrix(maps_out$unit.classif) 


##-------------------Shouldn't need to change much below this line --------------------
#pre-allocate matrix of m years (rows) and n nodes (columns)
maps_z <- matrix(data=NA,nrow=length(maps_uclass),ncol=nodes)
for (entry in c(1:length(maps_uclass))){
  maps_z[entry,maps_uclass[entry]]<-1
}

#generate Hovmoller plot
dev.new()
image(c(1:nodes), c(yearmin:yearmax), t(maps_z), 
      ylim = c(yearmax, yearmin),
      col = "black", ylab = "Year", xlab = "Node", main = title_name)

