##runSOM_xyf_super.r
#Emily Norton, 8/7/18
#This script uses the som, xyf, or (coming soon!) supersom function of the kohonen package to identify self-organizing maps
#for one, two, or three (or more!) variables at a time
#This program can be used to systematically run and test various SOMs, 
#and then generate several figures and/or map outputs based on the codes produced

#clear any lingering objects
rm(list = ls(all = TRUE))
#clear any potentially interfering loaded packages
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)

#Set the following options before running
setwd('/Users/emilynorton/Documents/FATE_Hunsicker_Bond_OcnAtmCoupling/RScripts/')

#option to set seed for reproducible results
set.seed(10)   

# Set number of variables we're going to compare 
compsize <- 1   #1 = univariable "som", 2 = two co-variable "xyf", 3 or more = multivariable "supersom"

#set variable filenames to load - only has data for the sea points
#NOTE: these files must contain headers and 'year' column, and must *at least* contain the temporal range of interest
#Non-spatially explicit data:
#file1 <- '../BiologyData/goa.biology.pca.1965.2012.csv'        # Contains V1 and V2
#file1 <- '../BiologyData/goa.biology.pca.1965.2012_V1only.csv'  # Contains V1 only
#file1 <- '../BiologyData/groundfish1.csv'
#file2 <- '../BiologyData/groundfish3.csv'
#file1 <- '../PhysicalData/pdo_1900to2018_FebtoApr_avg.csv'     # PDO averaged Feb to Apr
#file4 <- '../PhysicalData/npgo_1950to2018_FebtoApr_avg.csv'     # NPGO averaged Feb to Apr

#re-shaped spatially explicit data
file1 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds20to65N_120to255E_MLD.csv'   #MLD only available 1959-2011
#file4 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds40to65N_170to235E_MLD.csv'  
#file2 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds20to65N_120to255E_skt_detrend_incLat.csv'  #skt avail 1948 - 2018
#file1 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds40to65N_170to235E_skt_detrend_incLat.csv'
#file1 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds20to65N_120to255E_slp_incLat.csv'  #slp avail 1948 - 2018
#file1 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds40to65N_170to235E_slp_incLat.csv'
#file2 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds20to65N_120to255E_UFLX_Negated.csv'    #UFLX avail 1948 - 2018
#file2 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds40to65N_170to235E_UFLX_Negated.csv' 
#file3 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds20to65N_120to255E_VFLX_Negated.csv'    #VFLX avail 1948 - 2018
#file3 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds40to65N_170to235E_VFLX_Negated.csv'  
#file2 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds20to65N_120to255E_CURL_Negated.csv'
#file2 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds40to65N_170to235E_CURL_Negated.csv'
#file2 <- '../PhysicalData/Averages_Reshaped/MonthlyAvg11to3_forYears1948to2018_Bounds20to65N_130to255E_CURL_Negated.csv'


# Set variable grid sizes: "SLP"/"HGT" (both same grid), "SKT", "MLD", "UFLX", "VFLX", "CURL". If strictly timeseries var (e.g. biol V1), which doesn't have
# a spatial grid, it doesn't matter what you choose - the slp grid will be used by default - so just use the variable name
vargrid1 <- "MLD"
vargrid2 <- "SLP" 
vargrid3 <- "VFLX"
vargrid4 <- "MLD"

# Set temporal range of interest  **THIS IS NEW!! And allows you to re-use input files with longer temporal ranges
yearmin <- 1959
yearmax <- 2011

# Set spatial boundaries for the map grid  (seasonal averaged files must already be generated with these bounds - use matlab script "Average_PickYourMonths_SpatialBounds.m"
# and then add a "year" column and header row)
latmin = 40
latmax = 65
lonmin = 170
lonmax = 235

#set dimensions for the som grid (2D) and shape of somgrid (options: 'hexagonal','rectangular')
sdim1 <- 2 
sdim2 <- 1
sshape <- 'hexagonal' 

# set relative weights for each variable (up to 4 variables)
wei1 <- 1   #weight for var 1
wei2 <- 1   #weight for var 2
wei3 <- 1   #weight for var 3
wei4 <- 1   #weight for var 4


# Which types of plots do you want to make? 
booCodes <- 'T'
booQuality <- 'F'
booCounts <- 'F'
booYearNode <- 'T'
booGeoMap <- 'T'

# Plot title name and name to save plots
title_name <- ""   
plot_fname <- "test"

#Do you want to save SOM codes or maps as .csv files?
saveCodes <- 'F'
saveMaps <- 'F'

# Code and map filenames  ###COME BACK TO THIS>>>need to save depending on number of vars, etc.
codes1_fname <- sprintf("codes_SOM_%iby%i%sgrid_%ito%iyears_%s.csv", sdim1, sdim2, sshape, yearmin, yearmax, vargrid1)
codes2_fname <- sprintf("codes_SOM_%iby%i%sgrid_%ito%iyears_%s.csv", sdim1, sdim2, sshape, yearmin, yearmax, vargrid2)
codes3_fname <- sprintf("codes_SOM_%iby%i%sgrid_%ito%iyears_%s.csv", sdim1, sdim2, sshape, yearmin, yearmax, vargrid3)
codes4_fname <- sprintf("codes_SOM_%iby%i%sgrid_%ito%iyears_%s.csv", sdim1, sdim2, sshape, yearmin, yearmax, vargrid4)

varnames <- paste(c(vargrid1, vargrid2, vargrid3, vargrid4),collapse="_")
map_fname = sprintf("maps_unitclassif_SOM_%iby%i%sgrid_%ito%iyears_%s.csv", sdim1, sdim2, sshape, yearmin, yearmax, varnames)


##----------Shouldn't need to change too much below this line ----------------
nodes <- sdim1*sdim2  #this is the number of nodes that all codes are grouped into

#Load kohonen package (load other packages later)
library(tidyverse)   #for dplyr and ggplot - load first because of conflicts
library(kohonen)  #for SOMs


#Load csv file(s) and grab only the values of interest for SOM comparison (i.e. not headers)

data1 <- read_csv(file1)
data1_forSOM <- data1 %>%   #get rid of year column...may want to prepare in other ways, e.g. make sure all of the same years are used for both files, etc.
  filter(year <= yearmax, year >= yearmin) %>%
  select(-year)


if (compsize > 1) {
  data2 <- read_csv(file2)
  data2_forSOM <- data2 %>% 
    filter(year <= yearmax, year >= yearmin) %>%
    select(-year)
} 

if (compsize > 2){
  data3 <- read_csv(file3)
  data3_forSOM <- data3 %>% 
    filter(year <= yearmax, year >= yearmin) %>%
    select(-year)
}

if (compsize > 3){
  data4 <- read_csv(file4)
  data4_forSOM <- data4 %>% 
    filter(year <= yearmax, year >= yearmin) %>%
    select(-year)
}


#set grid filenames to load, with only lat/lon for the sea points, and the indexing info for each point - Updated 12/12/10 - reshaped grid
# May need to change these file paths to match your directory structure
#slpseaindsfile <- '../PhysicalData/slp_seainds.csv' 
#slplatvecfile <- '../PhysicalData/slp_lat_vec.csv'
#slplonvecfile <- '../PhysicalData/slp_lon_vec.csv'

#sktseaindsfile <- '../PhysicalData/skt_seainds.csv'
#sktlatvecfile <- '../PhysicalData/skt_lat_vec.csv'
#sktlonvecfile <- '../PhysicalData/skt_lon_vec.csv'

sktseaindsfile <- sprintf('../PhysicalData/Averages_Reshaped/skt_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE_incLat.csv', latmin, latmax, lonmin, lonmax) 
sktlatvecfile <- '../PhysicalData/Averages_Reshaped/skt_lat_vec_Bounds_OrigGrid_20to65N_120to255E_incLat.csv'
sktlonvecfile <- '../PhysicalData/Averages_Reshaped/skt_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv'

slpseaindsfile <- sprintf('../PhysicalData/Averages_Reshaped/slp_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE_incLat.csv', latmin, latmax, lonmin, lonmax) 
slplatvecfile <- '../PhysicalData/Averages_Reshaped/slp_lat_vec_Bounds_OrigGrid_20to65N_120to255E_incLat.csv'
slplonvecfile <- '../PhysicalData/Averages_Reshaped/slp_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv'

uflxseaindsfile <- sprintf('../PhysicalData/Averages_Reshaped/UFLX_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE.csv', latmin, latmax, lonmin, lonmax) 
uflxlatvecfile <- '../PhysicalData/Averages_Reshaped/UFLX_lat_vec_Bounds_OrigGrid_20to65N_120to255E.csv'
uflxlonvecfile <- '../PhysicalData/Averages_Reshaped/UFLX_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv'

vflxseaindsfile <- sprintf('../PhysicalData/Averages_Reshaped/VFLX_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE.csv', latmin, latmax, lonmin, lonmax) 
vflxlatvecfile <- '../PhysicalData/Averages_Reshaped/VFLX_lat_vec_Bounds_OrigGrid_20to65N_120to255E.csv'
vflxlonvecfile <- '../PhysicalData/Averages_Reshaped/VFLX_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv'

curlseaindsfile <- sprintf('../PhysicalData/Averages_Reshaped/CURL_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE.csv', latmin, latmax, lonmin, lonmax) 
curllatvecfile <- '../PhysicalData/Averages_Reshaped/CURL_lat_vec_Bounds_OrigGrid_20to65N_120to255E.csv'
curllonvecfile <- '../PhysicalData/Averages_Reshaped/CURL_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv'


mldseaindsfile <- sprintf('../PhysicalData/Averages_Reshaped/MLD_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE.csv', latmin, latmax, lonmin, lonmax) 
mldlatvecfile <- '../PhysicalData/Averages_Reshaped/MLD_lat_vec_Bounds_OrigGrid_20to65N_120to255E.csv'
mldlonvecfile <- '../PhysicalData/Averages_Reshaped/MLD_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv'

glist <- c(vargrid1, vargrid2, vargrid3, vargrid4)
grids <- glist[1:compsize]

#if (sum(grepl("SLP", grids, fixed=TRUE))>0){                     #always load these, since they are the default for non-spatial data
slpseainds <- as.matrix(read.csv(slpseaindsfile, header = FALSE)) 
slplatvec <-as.matrix(read.csv(slplatvecfile, header=FALSE))
slplonvec <-as.matrix(read.csv(slplonvecfile, header=FALSE))
#}

if (sum(grepl("SKT", grids, fixed=TRUE))>0){
sktseainds <- as.matrix(read.csv(sktseaindsfile, header = FALSE)) 
sktlatvec <-as.matrix(read.csv(sktlatvecfile, header=FALSE))
sktlonvec <-as.matrix(read.csv(sktlonvecfile, header=FALSE))
}

if (sum(grepl("UFLX", grids, fixed=TRUE))>0){
uflxseainds <- as.matrix(read.csv(uflxseaindsfile, header = FALSE)) 
uflxlatvec <-as.matrix(read.csv(uflxlatvecfile, header=FALSE))
uflxlonvec <-as.matrix(read.csv(uflxlonvecfile, header=FALSE))
}

if (sum(grepl("VFLX", grids, fixed=TRUE))>0){
vflxseainds <- as.matrix(read.csv(vflxseaindsfile, header = FALSE)) 
vflxlatvec <-as.matrix(read.csv(vflxlatvecfile, header=FALSE))
vflxlonvec <-as.matrix(read.csv(vflxlonvecfile, header=FALSE))
}

if (sum(grepl("CURL", grids, fixed=TRUE))>0){
  curlseainds <- as.matrix(read.csv(curlseaindsfile, header = FALSE)) 
  curllatvec <-as.matrix(read.csv(curllatvecfile, header=FALSE))
  curllonvec <-as.matrix(read.csv(curllonvecfile, header=FALSE))
}

if (sum(grepl("MLD", grids, fixed=TRUE))>0){
mldseainds <- as.matrix(read.csv(mldseaindsfile, header = FALSE)) 
mldlatvec <-as.matrix(read.csv(mldlatvecfile, header=FALSE))
mldlonvec <-as.matrix(read.csv(mldlonvecfile, header=FALSE))
}


# Load grid files in as matrices, and the indices ('seainds') where the data are from - NOTE: 'SLP' (sea level pressure)
# and 'HGT' (geopotential height at 200mbar) have the same grid and sea indices
if (vargrid1 == 'SKT') {
seaindsM1 <- sktseainds
latvecM1 <- sktlatvec
lonvecM1 <- sktlonvec
} else if (vargrid1 == 'UFLX') {
  seaindsM1 <- uflxseainds
  latvecM1 <- uflxlatvec
  lonvecM1 <- uflxlonvec
} else if (vargrid1 == 'VFLX') {
  seaindsM1 <- vflxseainds
  latvecM1 <- vflxlatvec
  lonvecM1 <- vflxlonvec
} else if (vargrid1 == 'CURL') {
  seaindsM1 <- curlseainds
  latvecM1 <- curllatvec
  lonvecM1 <- curllonvec
} else if (vargrid1 == 'MLD') {
  seaindsM1 <- mldseainds
  latvecM1 <- mldlatvec
  lonvecM1 <- mldlonvec
} else {                       # if vargrid 1 == "SLP" | vargrid2 == 'HGT' (or anything else, assume slp-sized grid)
seaindsM1 <- slpseainds
latvecM1 <- slplatvec
lonvecM1 <- slplonvec
}

if (compsize > 1) {
  if (vargrid2 == 'SKT') {
    seaindsM2 <- sktseainds
    latvecM2 <- sktlatvec
    lonvecM2 <- sktlonvec
  } else if (vargrid2 == 'UFLX') {
    seaindsM2 <- uflxseainds
    latvecM2 <- uflxlatvec
    lonvecM2 <- uflxlonvec
  } else if (vargrid2 == 'VFLX') {
    seaindsM2 <- vflxseainds
    latvecM2 <- vflxlatvec
    lonvecM2 <- vflxlonvec
  } else if (vargrid2 == 'CURL') {
    seaindsM2 <- curlseainds
    latvecM2 <- curllatvec
    lonvecM2 <- curllonvec
  } else if (vargrid2 == 'MLD') {
    seaindsM2 <- mldseainds
    latvecM2 <- mldlatvec
    lonvecM2 <- mldlonvec
  } else {
    seaindsM2 <- slpseainds
    latvecM2 <- slplatvec
    lonvecM2 <- slplonvec
  }
}

if (compsize > 2) {
  if (vargrid3 == 'SKT') {
    seaindsM3 <- sktseainds
    latvecM3 <- sktlatvec
    lonvecM3 <- sktlonvec
  } else if (vargrid3 == 'UFLX') {
    seaindsM3 <- uflxseainds
    latvecM3 <- uflxlatvec
    lonvecM3 <- uflxlonvec
  } else if (vargrid3 == 'VFLX') {
    seaindsM3 <- vflxseainds
    latvecM3 <- vflxlatvec
    lonvecM3 <- vflxlonvec
  } else if (vargrid3 == 'CURL') {
    seaindsM3 <- curlseainds
    latvecM3 <- curllatvec
    lonvecM3 <- curllonvec
  } else if (vargrid3 == 'MLD') {
    seaindsM3 <- mldseainds
    latvecM3 <- mldlatvec
    lonvecM3 <- mldlonvec
  } else {
    seaindsM3 <- slpseainds
    latvecM3 <- slplatvec
    lonvecM3 <- slplonvec
  }
}

if (compsize > 3) {
  if (vargrid4 == 'SKT') {
    seaindsM4 <- sktseainds
    latvecM4 <- sktlatvec
    lonvecM4 <- sktlonvec
  } else if (vargrid4 == 'UFLX') {
    seaindsM4 <- uflxseainds
    latvecM4 <- uflxlatvec
    lonvecM4 <- uflxlonvec
  } else if (vargrid4 == 'VFLX') {
    seaindsM4 <- vflxseainds
    latvecM4 <- vflxlatvec
    lonvecM4 <- vflxlonvec
  } else if (vargrid4 == 'CURL') {
    seaindsM4 <- curlseainds
    latvecM4 <- curllatvec
    lonvecM4 <- curllonvec
  } else if (vargrid4 == 'MLD') {
    seaindsM4 <- mldseainds
    latvecM4 <- mldlatvec
    lonvecM4 <- mldlonvec
  } else {
    seaindsM4 <- slpseainds
    latvecM4 <- slplatvec
    lonvecM4 <- slplonvec
  }
}

#Run SOMs with 1 ("som"), 2 ("xyf"), or 3 or more variables ("supersom") taken into account together
if (compsize == 1) {
  som_out <- som(scale(data1_forSOM),grid=somgrid(sdim1,sdim2,sshape))
}

if (compsize == 2) {
  som_out <- xyf(X=scale(data1_forSOM),Y=scale(data2_forSOM),user.weights=c(wei1,wei2),grid=somgrid(sdim1,sdim2,sshape))
}

if (compsize == 3) {
  Comb_dat <- list(scale(data1_forSOM), scale(data2_forSOM), scale(data3_forSOM))
  som_out <- supersom(Comb_dat, user.weights=c(wei1,wei2,wei3),grid=somgrid(sdim1,sdim2,sshape))
}

if (compsize == 4) {
  Comb_dat <- list(scale(data1_forSOM), scale(data2_forSOM), scale(data3_forSOM), scale(data4_forSOM))
  som_out <- supersom(Comb_dat, user.weights=c(wei1,wei2,wei3,wei4),grid=somgrid(sdim1,sdim2,sshape))
}

#Get codes for the SOM
codes_out <- getCodes(som_out)

#Get maps for the SOM
maps_out <- map(som_out)


# Save codes and maps, if desired
if (saveCodes == 'T') {
  write.table(codes_out1, file = codes1_fname, row.names=FALSE,col.names=FALSE,sep=',')
    if (compsize > 1) {
      write.table(codes_out2, file = codes2_fname, row.names=FALSE,col.names=FALSE,sep=',') }
    if (compsize > 2) {
      write.table(codes_out3, file = codes3_fname, row.names=FALSE,col.names=FALSE,sep=',') }
    if (compsize > 3) {
      write.table(codes_out4, file = codes4_fname, row.names=FALSE,col.names=FALSE,sep=',') }
}

if (saveMaps == 'T') {
write.table(maps_out$unit.classif, file = map_fname, row.names=FALSE,col.names=FALSE,sep=',')
}


#### Now begin making figures --------------------------------------
##Load libraries for making figures
#library(plotly)   #for plotting maps
library(maps)
library(mapproj)


#Make hovmoller-esque plot for what years were sorted into each node
if (booYearNode == 'T') {
    source('plot_year_node.r')
}

#Generate informational plots for kohonen
  source('plot_kohonen_fig.r')


#Generate geographic map of averaged codes for each node
if (booGeoMap == 'T') {
    source('plot_geo_map.r')
}
