library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(strucchange)

# data & code for reproducing Fig. 1, Changes in North Pacific atmosphere and ocean climate after 1988/89.

# start by loading NE Pacific SST
nc <- nc_open("data/fig1.sst.nc")

# get lat/long
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# assign dates
raw <- ncvar_get(nc, "time") # seconds since January 1, 1970
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# save month and year
m <- months(d)
yr <- years(d)

# get required sst data
SST <- ncvar_get(nc, "sst")

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix

# plot average temp at each cell to check
z <- colMeans(SST)
z <- t(matrix(z, length(y)))
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# looks good

# set names for matrix
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

f <- function(x) tapply(x, m[yr %in% 1951:1980], mean)  # function to compute monthly means for a single time series

mu <- apply(SST[yr %in% 1951:1980,], 2, f)	# Compute monthly means for each cell

# change to matrix
mu <- mu[rep(1:12, floor(length(d)/12)),] 

# add trailing months from incomplete 2019 data
xtra <- 12*((length(d)/12)-floor(length(d)/12))

mu <- rbind(mu, mu[1:xtra,])

# Now weight by cell area
weights <-  sqrt(cos(lat.t*pi/180))
ff <- function(x) weighted.mean(x, w=weights, na.rm=T)

# Calculate anomalies
anom <- SST-mu

# And get weighted mean anomaly for each month 
anom <- apply(anom,1,ff)   

# plot against decimal year
dec.yr <- as.numeric(as.character(yr)) + (as.numeric(m)-0.5)/12

plot(dec.yr, anom, type="l")
abline(h=0)

# now get annual winter (Nov-Mar) means, with year corresponding to January
# first, set winter year
win.yr <- as.numeric(as.character(yr))
# and advance for Nov / Dec
win.yr[m %in% c("Nov", "Dec")] <- win.yr[m %in% c("Nov", "Dec")]+1

# select anomalies only for these months
win.anom <- anom[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
win.yr <- win.yr[m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

win.anom <- tapply(win.anom, win.yr, mean)

# and plot to check
plot(1951:2019, win.anom[2:70], type="l", xlab="", ylab="Anomaly wrt 1951-1980, ºC")
abline(h=0)

# set up as time series for breakpoint test
sst.dat <- ts(data=win.anom[2:70], 1951, 2019, frequency=1)

# fit breakpoint model
bp.sst <- breakpoints(sst.dat ~ 1)
summary(bp.sst) 
# best model (lowest BIC) identifies 2 breaks
# first break is in 1998, second break in 2008 
# there is no obvious break in the time series in 2008 - I think we don't have enough data at the end of the time series
# to fit the break at the end of the time series, so I will reject the 2008 break and fit a model to a change at 2014
# note that this is a subjective decision!

# now fit a shifting mean model
# create a data frame
plot.dat <- data.frame(year=1951:2019, anom=win.anom[2:70])
# define eras for model: 1951-1988, 1989-2013, 2014-2019 (the latter is set to correspond to the onset of heatwaves in the N. Pacific)
plot.dat$era <- ifelse(plot.dat$year <= 1988, 1, ifelse(plot.dat$year %in% 1989:2013, 2, 3))

# fit a linear model to the eras
mod <- lm(anom ~ era, data=plot.dat)

# and get predicted values of the model to plot
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit

# now save for a combined plot

# set the colors to use - colorblind pallette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_classic())

a.plot <- ggplot(plot.dat, aes(year, anom)) +
  geom_line(size=0.4) +
  geom_line(aes(year, mean), color=cb[6], size=0.8) + 
  theme(axis.title.x = element_blank()) + ylab("Anomaly (ºC)") + ggtitle("Winter SST") +
  geom_vline(xintercept = 1988.5, lty=2) +
  xlim(1950,2020)

##################
# now SLP variability in the center of variability for the Aleutian Low

# first, load data
nc.slp <- nc_open("data/fig1.slp.nc")

# process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "time")  # seconds since 1-1-1970
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract area of interest
# 54-62 deg. N, 200-226 deg. E
x <- ncvar_get(nc.slp, "longitude")
y <- ncvar_get(nc.slp, "latitude")

# extract data
SLP <- ncvar_get(nc.slp, "slp", verbose = F)

# Change data to a matrix
SLP <- aperm(SLP, 3:1)  

# Reverse order of latitudes to be increasing for plotting
SLP <- SLP[,5:1,]  

# Also reverse  vector of latitudes
y <- rev(y)  

# Change to matrix
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Get lat/long vectors and add names to SLP matrix
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
z <- colMeans(SLP)
z <- t(matrix(z, length(y)))
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", ylim=c(35,66), xlim=c(170,220))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# looks good

# now we need to get monthly means!
m.y <- paste(years(d), as.numeric(months(d)), sep="-") # make a month-year factor from the dates

f <- function(x) tapply(x, m.y, mean)
SLP.m <- as.data.frame(apply(SLP, 2, f))

vv <- matrix(unlist(strsplit(as.character(rownames(SLP.m)), "-")),ncol=2, byrow = T)

SLP.m$year <- as.numeric(vv[,1])
SLP.m$month <- as.numeric(vv[,2])

SLP.m <- SLP.m %>%
  arrange(month) %>%
  arrange(year)

# remove seasonal signal
m <- SLP.m$month
yr <- SLP.m$year

f <- function(x) tapply(x, m, mean)
mu <- apply(SLP.m, 2, f)	# Compute monthly means for each cell

# process as for SST
mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

SLP.anom <- SLP.m[,1:35] - mu   # Compute matrix of anomalies - dropping year and month!

# get average anomaly across the area
SLP.anom <- rowMeans(SLP.anom)

# fit to winter means
win.yr <- ifelse(m %in% c(11,12), yr+1, yr)
SLP.win.anom <- SLP.anom[m %in% c(11,12,1:3)]
win.yr <- win.yr[m %in% c(11,12,1:3)]
SLP.win.anom <- tapply(SLP.win.anom, win.yr, mean)
plot(1949:2019, SLP.win.anom, type="l")

# and calculate standard deviation over 11-month rolling windows
SLP.win.sd <- rollapply(SLP.win.anom, 11, sd, fill=NA)
plot(1949:2019, SLP.win.sd, type="l")

# now fit a non-parametric regression

# first, make a data frame
plot.dat <- data.frame(year=1954:2014, sd=na.omit(SLP.win.sd))

# fit the model
mod <- gam(sd ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit                    

b.plot <- ggplot(plot.dat, aes(year, sd)) +
  geom_line(size=0.4) +
  geom_line(aes(year, mean), color=cb[6], size=0.8) + # geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.2) +
  theme(axis.title.x = element_blank()) + ylab("Standard deviation (pa)") + ggtitle("Aleutian Low variability") +
  geom_vline(xintercept = 1988.5, lty=2) +
  xlim(1950,2020)

##########################
# finally, pdo-npgo correlation
library(tidyr)

# load pdo
# download.file("http://jisao.washington.edu/pdo/PDO.latest", "data/pdo") # uncomment to load
names <- read.table("data/pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("data/pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>% # not loading tidyr because I think it conflicts with maps!
  arrange(YEAR)

# load npgo
# download.file("http://www.oces.us/npgo/npgo.php", "data/npgo") # uncomment to load
npgo <- read.table("data/npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# set up data frame and loop through the rolling windows
pdo.npgo <- data.frame(year=rep(1950:2018, each=12), month=rep(1:12, 69), npgo=npgo$value, pdo=pdo$value[601:1428], cor=NA)

for(i in 67:(nrow(pdo.npgo)-66)){ # 132 month rolling windows!
  
  pdo.npgo$cor[i] <- cor(pdo.npgo$npgo[(i-66):(i+66)], pdo.npgo$pdo[(i-66):(i+66)])
}

# add decimal year for plotting
pdo.npgo$dec.yr <- pdo.npgo$year + (pdo.npgo$month-0.5)/12

# fit gam - limiting effective degrees of freedom with k=8
mod <- gam(cor ~ s(dec.yr, k=8), data=pdo.npgo, correlation=corAR1())
pred <- predict(mod, se=T)

pdo.npgo$fit <- c(rep(NA, 66), pred$fit,rep(NA, 828-759))

c.plot <- ggplot(pdo.npgo, aes(dec.yr, cor)) +
  theme_classic() +
  # geom_hline(yintercept = 0, color="dark grey") +
  geom_line(size=0.4) +
  geom_line(aes(dec.yr, fit), color=cb[6], size=0.8) + # geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.2) +
  xlab("") + ylab("Pearson's correlation") + ggtitle("PDO-NPGO correlation") +
  xlim(1950,2020) +
  geom_vline(xintercept = 1988.5, lty=2)

png("Fig 1.png", 4, 7, units="in", res=300) 
ggarrange(a.plot, b.plot, c.plot, labels = c("a)", "b)", "c)"),  nrow=3, align="v")
dev.off()








######
# first, load all the data - SLP,SST, PDO, NPGO

# load monthly NCEP/NCAR SLP data
nc.slp <- nc_open("data/fig1.slp.nc")

# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "TIME")  # seconds since 1-1-1970
d <- dates(raw, origin = c(1,1,0001))
yr <- years(d)

# define lat/long
x <- ncvar_get(nc.slp, "LON53_101")
y <- ncvar_get(nc.slp, "LAT45_69")

# and extract data
SLP <- ncvar_get(nc.slp, "SLP", verbose = F)

# Change to a matrix for ease of handling
# Reverse order
SLP <- aperm(SLP, 3:1)  

# And change into a matrix (rows = months, columns = cells)
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Get vectors of latitudes and longitudes
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   

# Add dimnames for the matrix
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
z <- colMeans(SLP) # Mean cell values as data to plot
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F,add=T, lwd=1)
#looks good!

# load pdo
# download.file("http://jisao.washington.edu/pdo/PDO.latest", "data/pdo") # uncomment to load
names <- read.table("data/pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("data/pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  tidyr::gather(month, value, -YEAR) %>% # not loading tidyr because I think it conflicts with maps!
  arrange(YEAR)

# load npgo
# download.file("http://www.oces.us/npgo/npgo.php", "data/npgo") # uncomment to load
npgo <- read.table("data/npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# load ERSSTv5 SST data
# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2019-03-01)][(0.0):1:(0.0)][(30):1:(66)][(150):1:(250)]", "data/fig1.sst.nc")
# uncomment that line to download the data

nc <- nc_open("data/fig1.sst.nc")

# get lat/long
x.t <- ncvar_get(nc, "longitude")
y.t <- ncvar_get(nc, "latitude")
lat.t <- rep(y.t, length(x.t))   # Vector of latitudes
lon.t <- rep(x.t, each = length(y.t))   # Vector of longitudes

# assign dates
raw <- ncvar_get(nc, "time") # seconds since January 1, 1970
h <- raw/(24*60*60)
d.t <- dates(h, origin = c(1,1,1970))

# year for processing later
m <- as.numeric(months(d.t))
yr <- as.numeric(as.character(years(d.t)))

# get required sst data
SST <- ncvar_get(nc, "sst")

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix

# plot to check
z <- colMeans(SST)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# looks good

# set names
dimnames(SST) <- list(as.character(d.t), paste("N", lat.t, "E", lon.t, sep=""))

###########################
# calculate the various time series
###########################

# First, mean winter temp

# Create a function to compute monthly means for a single time series for the 1951-1980 climatology
f <- function(x) tapply(x, m[yr %in% 1951:1980], mean)  

# Compute monthly means for each cell
mu <- apply(SST[yr %in% 1951:1980,], 2, f)	

# Change into a matrix with 12 monthly values per year
mu <- mu[rep(1:12, floor(length(d.t)/12)),] 

# Add the trailing months that come from the final incomplete year
xtra <- 12*((length(d.t)/12)-floor(length(d.t)/12))
mu <- rbind(mu, mu[1:xtra,])

# Now weight by cell area
weights <-  sqrt(cos(lat.t*pi/180))
ff <- function(x) weighted.mean(x, w=weights, na.rm=T)

# Calculate anomalies
SST.anom <- SST-mu

# And get weighted mean anomaly for each month 
SST.anomTS <- apply(SST.anom,1,ff)   

# Set decimal year for plotting
dec.yr.t <- as.numeric(as.character(yr)) + (as.numeric(m)-0.5)/12

# Plot to check
plot(dec.yr.t, SST.anomTS, type="l")
abline(h=0)

# Now, calculate time series of SLP anomalies for Aleutian Low area

# First, subset SLP to AL area...45-55 N, 192.5-207.5 E
latAL <- lat >= 46.25 & lat <=56.25
lonAL <- lon >= 186.25 & lon <=206

SLP.AL <- SLP
SLP.AL[,!latAL] <- NA
SLP.AL[,!lonAL] <- NA

# Plot to check
z <- colMeans(SLP.AL)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
# Looks good!

# Calculate monthly anomalies to remove seasonal signal - same approach as used for SST
m <- months(d)
yr <- years(d)

f <- function(x) tapply(x, m, mean)
mu <- apply(SLP.AL, 2, f)	

mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

SLP.ALanom <- SLP.AL - mu   # Compute matrix of anomalies - dropping year and month!

# get average anomaly across the area
SLP.ALanom <- rowMeans(SLP.ALanom, na.rm=T)

# smooth with 11-m rolling mean
SLP.sm <- rollmean(SLP.ALanom,11,fill=NA)

# and get the rolling 21-yr (253-mo) sd
SLP.sd <- rollapply(SLP.sm, 253, sd, fill=NA)

# Plot against decimal year to check
dec.yr <- as.numeric(as.character(yr))+(as.numeric(m)-0.5)/12
plot(dec.yr, SLP.sm, type="l")

# now pdo-npgo correlations for Jan. 1950 - Sept. 2018
pdoTS <-pdo[601:1425,]
npgoTS <- npgo[1:825,]
rownames(pdoTS) <- rownames(npgoTS) <- 1:nrow(pdoTS)

# make an object to catch the rolling wondow correlations
pdo.npgo <- NA

for(i in 127:(nrow(pdoTS)-126)){
  # i <- nrow(pdoTS)-126
  pdo.npgo[(i-126)] <- cor(pdoTS$value[(i-126):(i+126)], npgoTS$value[(i-126):(i+126)])
  
}










# put together time series for DFA and plots

dfa.dat <- data.frame(dec.yr=dec.yr, year=as.numeric(as.character(yr)), month=as.numeric(m), sst=SST.anom[match(dec.yr, dec.yr.t)], AL.sd=SLP.sd, PDO.NPGO.cor=NA, SLP.PDO.NS=NA, SLP.NPGO=NA)
plot(dfa.dat$dec.yr,dfa.dat$sst, type="l")
plot(dec.yr.t, SST.anom, type="l")

# and get full-field SD values by era to plot

# recalculate anomalies for the entire era
mu <- apply(SLP, 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

SLPanom <- SLP - mu   # Compute matrix of anomalies - dropping year and month!

ff <- function(x) rollmean(x,11,fill=NA)

# smooth with 11-m rolling mean
SLP.sm <- apply(SLPanom, 2, ff)

SLPsd1 <- apply(SLPanom[yr<=1988,], 2, sd, na.rm=T)
SLPsd2 <- apply(SLPanom[yr %in% 1989:2013,], 2, sd, na.rm=T)
SLPsd3 <- apply(SLPanom[yr >= 2014,], 2, sd, na.rm=T)

lim <- range(SLPsd1, SLPsd2)

xx <- c(186.25, 186.25, 206, 206, 186.25)
yy <- c(46.25, 56.25, 56.25, 46.25, 46.25)

# plot to check
 par(mfrow=c(1,2))


z <- SLPsd1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", zlim=lim)
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
mtext("1948-1988")
lines(xx,yy, lwd=1.5, col="magenta")

z <- SLPsd2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=lim)
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
mtext("1989-2013")
lines(xx,yy, lwd=1.5, col="magenta")


View(dfa.dat)
dfa.dat$PDO.NPGO.cor[151:723] <- pdo.npgo

# and attempt some plots that depict the declining independent predictive skill of the NPGO
# reload SST
nc <- nc_open("~updated.sst")

# get lat/long
x.t <- ncvar_get(nc, "longitude")
y.t <- ncvar_get(nc, "latitude")
lat.t <- rep(y.t, length(x.t))   # Vector of latitudes
lon.t <- rep(x.t, each = length(y.t))   # Vector of longitudes

# assign dates
raw <- ncvar_get(nc, "time") # seconds since January 1, 1970
h <- raw/(24*60*60)
d.t <- dates(h, origin = c(1,1,1970))

# year for processing later
m <- as.numeric(months(d.t))
yr <- as.numeric(as.character(years(d.t)))

# get required sst data
SST <- ncvar_get(nc, "sst")

# need to change to matrix for easier use
SST <- aperm(SST, 3:1) # transpose array

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  # Change to matrix

# plot to check
z <- colMeans(SST)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)

# set names
dimnames(SST) <- list(as.character(d.t), paste("N", lat.t, "E", lon.t, sep=""))


f <- function(x) tapply(x, m[yr %in% 1951:1980], mean)  # function to compute monthly means for a single time series

mu <- apply(SST[yr %in% 1951:1980,], 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d.t)/12)),] 

xtra <- 12*((length(d.t)/12)-floor(length(d.t)/12))

mu <- rbind(mu, mu[1:xtra,])

SST.anom <- SST-mu

# and remove land
land <- is.na(colMeans(SST.anom))
SST.anom <- SST.anom[,!land]

# separate into 1950:1988 and 1989:2013
SST.anom1 <- SST.anom[yr <= 1988,]
SST.anom2 <- SST.anom[yr %in% 1989:2013,]

# and separate pdo/npgo to same periods
pdo1 <- pdo$value[pdo$YEAR %in% 1950:1988]
pdo2 <- pdo$value[pdo$YEAR %in% 1989:2013]

npgo1 <- npgo$value[npgo$Year %in% 1950:1988]
npgo2 <- npgo$value[npgo$Year %in% 1989:2013]


# look at NPGO regression on SSTa
npgo.sst.regr1 <- npgo.sst.regr2 <- NA

for(i in 1:ncol(SST.anom1)){
  # i <- 1
  mod <- lm(SST.anom1[,i] ~ npgo1)
  npgo.sst.regr1[i] <- summary(mod)$coefficients[2,1]
  
  mod <- lm(SST.anom2[,i] ~ npgo2)
  npgo.sst.regr2[i] <- summary(mod)$coefficients[2,1]
}

# plot

lim <- range(npgo.sst.regr1, npgo.sst.regr2)

par(mfrow=c(2,2), mar=c(0.5, 0.5,2,2))
z <- rep(NA, ncol(SST))

z[!land] <- npgo.sst.regr1

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO regr. SSTa 1950-1988", cex=0.8)

z <- rep(NA, ncol(SST))
z[!land] <- npgo.sst.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO regr. EOF2 1989-2102", cex=0.8)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern1  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r1  # replace elements NOT corresponding to land with loadings!
z[!land] <- npgo.eof2.regr1  # replace elements NOT corresponding to land with loadings!

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO regr. EOF2 1950-1988", cex=0.8)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern2  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r2
z[!land] <- npgo.eof2.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO regr. EOF2 1989-2102", cex=0.8)

###
# and SLP-PDO
# reload SLP to make things easy...
nc.slp <- nc_open("/Users/MikeLitzow 1/Documents/R/climate scripts/64AA9803F7345DFE8991A731014FFB01_ferret_listing.nc")

# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "TIME")  # seconds since 1-1-1970
# h <- raw/(24*60*60)
d <- dates(raw, origin = c(1,1,0001))
yr <- years(d)


x <- ncvar_get(nc.slp, "LON53_101")
y <- ncvar_get(nc.slp, "LAT45_69")

# save to plot below
x.slp <- x
y.slp <- y

SLP <- ncvar_get(nc.slp, "SLP", verbose = F)
# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check
z <- colMeans(SLP)   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n")
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)

# load pdo
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load npgo
download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# smooth SLP with three month rolling mean
f <- function(x) rollmean(x,3,fill=NA)

SLPsm <- apply(SLP, 2, f)

rownames(SLPsm)

# limit to 1950:2017
SLPsm <- SLPsm[yr %in% 1950:2017,]

# limit PDO and NPGO to March 1950 - Feb 2018
# also expanding to split out 1950-1988 and 1989-2013
pdo[c(603,1418),]
pdo[c(603,1068),]
pdo[c(1069,1368),]

pdoTS <- pdo$value[603:1418]
pdoTS1 <- pdo$value[603:1068]
pdoTS2 <- pdo$value[1069:1368]

npgo[c(3,818),]
npgo[c(3,468),]
npgo[c(469,768),]

npgoTS <- npgo$value[3:818]
npgoTS1 <- npgo$value[3:468]
npgoTS2 <- npgo$value[469:768]

# and two era time series
rownames(SLPsm)[c(1,length(pdoTS1))]
rownames(SLPsm)[c(length(pdoTS1),length(pdoTS1)+length(pdoTS2))]
SLPsm1 <- SLPsm[1:length(pdoTS1),]
SLPsm2 <- SLPsm[(1+length(pdoTS1)):(length(pdoTS1)+length(pdoTS2)),]

# separate regressions in each era!
pdo.regr1 <- pdo.regr2 <- npgo.regr1 <- npgo.regr2 <- NA

for(i in 1:ncol(SLPsm)){
  #  i <- 1
  mod <- lm(SLPsm1[,i] ~ pdoTS1)
  pdo.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLPsm2[,i] ~ pdoTS2)
  pdo.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLPsm1[,i] ~ npgoTS1)
  npgo.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLPsm2[,i] ~ npgoTS2)
  npgo.regr2[i] <- summary(mod)$coef[2,1] 
}

# and plot
# set up color schemes
new.col <- my.col <- tim.colors(64)
grays <- c("gray98", "gray97", "gray96", "gray95", "gray94", "gray93", "gray92", "gray91", "gray90", "gray89", "gray88",
           "gray87", "gray86", "gray85", "gray84", "gray83", "gray82", "gray81", "gray80", "gray79", "gray78", "gray77",
           "gray76", "gray75", "gray74", "gray73", "gray72", "gray71", "gray71", "gray71", "gray70", "gray69", "gray68")


my.col[1:33] <- grays
my.col[22:43] <- c(grays[11:1], grays[1:11])

new.col[27:36] <- c(grays[5:1], grays[1:5])

lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)


# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(1.5,1.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,2), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("PDO forcing pattern 1950-1988", cex=0.8)

xx1 <- c(183.5, 183.5, 204, 204, 183.5)
yy1 <- c(46.5, 54, 54, 46.5, 46.5)

xx2 <- c(196, 196, 216.5, 216.5, 196)
yy2 <- c(39, 46.5, 46.5, 39, 39)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

# PDO second era
z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("PDO forcing pattern 1989-2013", cex=0.8)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")


# npgo first era
z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO forcing pattern 1950-1988", cex=0.8)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

# npgo second era
z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("NPGO forcing pattern 1989-2013", cex=0.8)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

###########################
# get rolling regression coefficients for each square

PDOkeep1 <- PDOkeep2 <- NPGOkeep <- NA
pdoSLP1 <- pdoSLP2 <- npgoSLP <- SLPsm

for(i in 1:length(lat)){
  # i <- 1
  PDOkeep1[i] <- inpolygon(lon[i], lat[i], xx1, yy1)
  PDOkeep2[i] <- inpolygon(lon[i], lat[i], xx2, yy2)
  NPGOkeep[i] <- inpolygon(lon[i], lat[i], xx3, yy3)
}

pdoSLP1[,!PDOkeep1] <- NA
pdoSLP2[,!PDOkeep2] <- NA
npgoSLP[,!NPGOkeep] <- NA


# plot to check
# first PDO1
z <- colMeans(pdoSLP1, na.rm=T)  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=tim.colors(64), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
lines(xx1, yy1, lwd=1.5, col="magenta")

# PDO2 
z <- colMeans(pdoSLP2, na.rm=T)  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=tim.colors(64), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
lines(xx2, yy2, lwd=1.5, col="magenta")

# now NPGO
z <- colMeans(npgoSLP, na.rm=T)  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image(x,y,z, col=tim.colors(64), ylim=c(20,68),
      xlab = "", ylab = "")

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
lines(xx3, yy3, lwd=1.5, col="magenta")

# all look good

# now try 21-yr rolling regressions!
pdo.regr1 <- pdo.regr2 <- npgo.regr <- NA

# get monthly means of each
pdoSLP1 <- rowMeans(pdoSLP1, na.rm=T)
pdoSLP2 <- rowMeans(pdoSLP2, na.rm=T)
npgoSLP <- rowMeans(npgoSLP, na.rm=T)

names(pdoSLP1)[1:(12*21)]
# will use 253 month (21 yr + 1 mo) rolling windows

roll <- data.frame(dec.yr=as.numeric(as.character(years(names(pdoSLP1))))+(as.numeric(months(names(pdoSLP1)))-0.5)/12, pdo.regr1=NA, pdo.regr2=NA, npgo.regr=NA)

# and the data for regressions
dat <- data.frame(year=as.numeric(as.character(years(names(pdoSLP1))))+(as.numeric(months(names(pdoSLP1)))-0.5)/12,
                  pdoSLP1=pdoSLP1, pdoSLP2=pdoSLP2, npgoSLP=npgoSLP, pdoTS=pdoTS, npgoTS=npgoTS, pdo.regr1=NA, pdo.regr2=NA, npgo.regr=NA) # so the correct lags are built in


for(i in 127:(nrow(dat)-126)){
  # i <- 127
  temp <- dat[(i-126):(i+126),]
  
  mod <- lm(temp$pdoSLP1 ~ temp$pdoTS)
  dat$pdo.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(temp$pdoSLP2 ~ temp$pdoTS)
  dat$pdo.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(temp$npgoSLP ~ temp$npgoTS)
  dat$npgo.regr[i] <- summary(mod)$coef[2,1]
  
}

dat$pdoNS <- dat$pdo.regr1/dat$pdo.regr2
  
dfa.dat$SLP.PDO.NS  <- dat$pdoNS[match(dfa.dat$dec.yr, dat$year)]
dfa.dat$SLP.NPGO  <- dat$npgo.regr[match(dfa.dat$dec.yr, dat$year)]

plot.dat <- dat %>%
  select(year, pdo.regr1, pdo.regr2, npgo.regr) %>%
  gather(key, value, -year)

ggplot(plot.dat, aes(year, value, color=key)) +
  theme_linedraw() +
  geom_line() +
  ylim(c(-0.5, -2.8)) + 
  xlim(c(1962,2005))

dat$pdo.ratio <- ifelse(is.na(dat$pdo.regr1),NA, dat$pdo.regr1/dat$pdo.regr2)

plot.dat <- dat %>%
  select(year, pdo.ratio, npgo.regr) %>%
  gather(key, value, -year)

ggplot(plot.dat, aes(year, value)) +
  theme_linedraw() +
  geom_line() +
  facet_wrap(~key, scales="free") + 
  xlim(c(1962,2005))

# now! just fit the DFA model...
all.dat <- as.matrix(t(dfa.dat[,4:8]))
colnames(all.dat) <- dfa.dat[,1]

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")

cntl.list = list(minit=200, maxit=10000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
model.data = data.frame()

# restrict to 1960-2008
some.dat <- all.dat[,colnames(all.dat) >=1960 & colnames(all.dat) <= 2008]

# # fit lots of models & store results
# # NOTE: this will take a long time to run!
# for(R in levels.R) {
#   for(m in 1:1) {  # allowing only 1 trend!
#     dfa.model = list(A="zero", R=R, m=m)
#     
# 
#     
#     kemz = MARSS(some.dat, model=dfa.model, control=cntl.list,
#                  form="dfa", z.score=TRUE)
#     model.data = rbind(model.data,
#                        data.frame(R=R,
#                                   m=m,
#                                   logLik=kemz$logLik,
#                                   K=kemz$num.params,
#                                   AICc=kemz$AICc,
#                                   stringsAsFactors=FALSE))
#     assign(paste("kemz", m, R, sep="."), kemz)
#   } # end m loop
# } # end R loop
# 
# model.data


# unconstrained hasn;t converged
# re-fit...
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)
model.list = list(A="zero", m=1, R="unconstrained") 
mod.check = MARSS(some.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)
# AICc is 4507.108... best!


# and plot
# get CI and plot loadings...
modCI <- MARSSparamCIs(mod.check)
modCI

new.names <- c("SST", "Aleutian Low SD", "PDO-NPGO correlation", "SLP-PDO North-South", "SLP-NPGO")
plot.CI <- data.frame(plot.names=new.names, names=rownames(all.dat), mean=modCI$par$Z, upCI=modCI$par.upCI$Z,
                           lowCI=modCI$par.lowCI$Z)

plot.CI <- arrange(plot.CI, mean)
plot.CI$names.order <- reorder(plot.CI$plot.names, plot.CI$mean)
dodge <- position_dodge(width=0.9)

png("basin climate DFA loadings.png", 4, 4, units="in", res=300)
ggplot(plot.CI, aes(x=names.order, y=mean)) + geom_bar(position="dodge", stat="identity") +
  theme_linedraw() + 
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  geom_hline(yintercept = 0) 
dev.off()

# save for combined plot
lp <- ggplot(plot.CI, aes(x=names.order, y=mean)) + geom_bar(position="dodge", stat="identity") +
  theme_linedraw() + 
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") + xlab("") + theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  geom_hline(yintercept = 0) + ggtitle("a) Loadings") 

# get CI for latent trends and plot
d <- tidy.marssMLE(mod.check, type="states")
# change time step to years
d$t <- as.numeric(colnames(some.dat))
# d$t <- 1:nrow(d)
# call this 'tp' (trend plot)
tp <- ggplot(data = d) +
  theme_linedraw() +
  geom_line(aes(t, estimate)) +  geom_hline(yintercept = 0) +
  theme(axis.title.x = element_blank()) + ylab("Trend") + ggtitle("b) Trend")


#####################
#####################

# sst time series plot
use.dat <- data.frame(yr=dec.yr.t, sst=SST.anomTS)
use.dat$era <- as.factor(ifelse(use.dat$yr<1989,1,2))

pred.sst1 <- predict(lm(sst ~ yr:era + era, data=use.dat))

pred.sst2 <- predict(gam(sst ~ s(yr), data=use.dat))

pred.sst3 <- rollmean(SST.anomTS, 133, fill=NA)

png("ne pacific sst plot.png", 5,4, units="in", res=300)
par(mar=c(4,4,1.5,1),las=1)
plot(dec.yr.t, SST.anomTS, type="l", xlab="", ylab="ºC wrt 1951-1980", xlim=xlim, col=cb[6])
lines(dec.yr.t, pred.sst2, col="red")
abline(h=0)
abline(v=1989, lty=2)
dev.off()



#######
# combined plot
png("reduced basin scale combined plot.png", 6,6, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my.col <- tim.colors(64)
grays <- c("gray98", "gray98", "gray97", "gray97", "gray96", "gray96", "gray95", "gray95", "gray94", "gray94", "gray93",
           "gray93", "gray92", "gray92", "gray91", "gray91", "gray90", "gray90", "gray89", "gray89", "gray88", "gray88",
           "gray87", "gray87", "gray86", "gray86", "gray85", "gray85", "gray84", "gray84", "gray83", "gray83", "gray82")


my.col[1:33] <- grays
# my.col[22:43] <- c(grays[11:1], grays)
# new.col[27:36] <- c(grays[5:1], grays[1:5])

xlim <- c(min(dec.yr), max(dec.yr.t))

par(mar=c(1.25,1.25,1.25,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(4,3), cex.axis=0.8, cex.lab=0.8, oma=c(0,2,2,0))


###
# now AL SD

lim <- range(SLPsd1, SLPsd2)

xx <- c(186.25, 186.25, 206, 206, 186.25)
yy <- c(46.25, 56.25, 56.25, 46.25, 46.25)


z <- SLPsd1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n",
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("a", adj=0.05, line=-1.4, cex=1)

z <- SLPsd2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n", 
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("b", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.5,3,1.5,1))

plot(dec.yr, SLP.sd, type="l", xlab="", ylab="Standard dev. (Pa)", xlim=xlim, col=cb[6])
abline(h=mean(SLP.sd, na.rm=T))
abline(v=1989, lty=2)
mtext("c", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.25,1.25,1.25,1))

# now atmospheric forcing of PDO/NPGO
lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)
# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("d", adj=0.05, line=-1.4, cex=1)

xx1 <- c(183.5, 183.5, 204, 204, 183.5)
yy1 <- c(46.5, 54, 54, 46.5, 46.5)

xx2 <- c(196, 196, 216.5, 216.5, 196)
yy2 <- c(39, 46.5, 46.5, 39, 39)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

# PDO second era
z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("e", adj=0.05, line=-1.4, cex=1)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$pdoNS, type="l", xlab="", ylab="North:South forcing ratio", xlim=xlim, col=cb[6])
abline(h=mean(dat$pdoNS, na.rm=T))
abline(v=1989, lty=2)
mtext("f", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.25,1.25,1.25,1))

# npgo first era
z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("g", adj=0.05, line=-1.4, cex=1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

# npgo second era
z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("h", adj=0.05, line=-1.4, cex=1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$npgo.regr, type="l", xlab="", ylab="Regression coef (Pa)", xlim=xlim, col=cb[6],
     ylim=c(max(dat$npgo.regr, na.rm=T), min(dat$npgo.regr, na.rm=T)))
abline(h=mean(dat$npgo.regr, na.rm=T))
abline(v=1989, lty=2)
mtext("i", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.25,1.25,1.25,1))


# finally, the NPGO-SST plots
lim <- range(npgo.sst.regr1, npgo.sst.regr2)
z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern1  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r1  # replace elements NOT corresponding to land with loadings!
z[!land] <- npgo.sst.regr1  # replace elements NOT corresponding to land with loadings!

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("j", adj=0.05, line=-1.4, cex=1)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern2  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r2
z[!land] <- npgo.sst.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("k", adj=0.05, line=-1.4, cex=1)

par(mar=c(1.5,3,1.5,1))

plot(dfa.dat$dec.yr, dfa.dat$PDO.NPGO.cor, type="l", xlab="", ylab="PDO-NPGO correlation", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("l", adj=0.05, line=-1.4, cex=1)
dev.off()


########
# older version
png("older larger basin scale combined plot.png", 7,10, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
 my.col <- tim.colors(64)
grays <- c("gray98", "gray98", "gray97", "gray97", "gray96", "gray96", "gray95", "gray95", "gray94", "gray94", "gray93",
           "gray93", "gray92", "gray92", "gray91", "gray91", "gray90", "gray90", "gray89", "gray89", "gray88", "gray88",
           "gray87", "gray87", "gray86", "gray86", "gray85", "gray85", "gray84", "gray84", "gray83", "gray83", "gray82")


my.col[1:33] <- grays
# my.col[22:43] <- c(grays[11:1], grays)
# new.col[27:36] <- c(grays[5:1], grays[1:5])

xlim <- c(min(dec.yr), max(dec.yr.t))

par(mar=c(1.5,1.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(5,3), cex.axis=0.8, cex.lab=0.8, oma=c(0,2,2,0))

lim <- range(SSTanom1, SSTanom2, na.rm=T)

z <- SSTanom1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]), 
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
mtext("a", adj=0.05, line=-1.4, cex=1.1)

z <- SSTanom2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=tim.colors(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]),
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x.t,y.t,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
mtext("b", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))
plot(dec.yr.t, SST.anomTS, type="l", xlab="", ylab="ºC wrt 1951-1980", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("c", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

###
# now AL SD

lim <- range(SLPsd1, SLPsd2)

xx <- c(186.25, 186.25, 206, 206, 186.25)
yy <- c(46.25, 56.25, 56.25, 46.25, 46.25)


z <- SLPsd1   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n",
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("d", adj=0.05, line=-1.4, cex=1.1)

z <- SLPsd2   # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=my.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=lim, ylim=c(20,70), yaxt="n", xaxt="n", 
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,80),add=T, lwd=1)
lines(xx,yy, lwd=2, col="magenta")
mtext("e", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))

plot(dec.yr, SLP.sd, type="l", xlab="", ylab="Standard dev. (Pa)", xlim=xlim, col=cb[6])
abline(h=mean(SLP.sd, na.rm=T))
abline(v=1989, lty=2)
mtext("f", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

# now atmospheric forcing of PDO/NPGO
lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)
# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
      xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("g", adj=0.05, line=-1.4, cex=1.1)

xx1 <- c(183.5, 183.5, 204, 204, 183.5)
yy1 <- c(46.5, 54, 54, 46.5, 46.5)

xx2 <- c(196, 196, 216.5, 216.5, 196)
yy2 <- c(39, 46.5, 46.5, 39, 39)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

# PDO second era
z <- pdo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
#image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("h", adj=0.05, line=-1.4, cex=1.1)

lines(xx1, yy1, lwd=1.5, col="magenta")
lines(xx2, yy2, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$pdoNS, type="l", xlab="", ylab="North:South forcing ratio", xlim=xlim, col=cb[6])
abline(h=mean(dat$pdoNS, na.rm=T))
abline(v=1989, lty=2)
mtext("i", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))

# npgo first era
z <- npgo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("j", adj=0.05, line=-1.4, cex=1.1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

# npgo second era
z <- npgo.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y.slp)))  # Convert vector to matrix and transpose for plotting
# image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
#            xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
image.plot(x.slp,y.slp,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68),
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x.slp,y.slp,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("k", adj=0.05, line=-1.4, cex=1.1)

xx3 <- c(194, 194, 199, 211, 216.5, 216.5, 194)
yy3 <- c(51.5, 54, 56.5, 61.5, 61.5, 51.5, 51.5)

lines(xx3, yy3, lwd=1.5, col="magenta")

par(mar=c(1.5,3,1.5,1))

plot(dat$year, dat$npgo.regr, type="l", xlab="", ylab="Regression coef (Pa)", xlim=xlim, col=cb[6],
     ylim=c(max(dat$npgo.regr, na.rm=T), min(dat$npgo.regr, na.rm=T)))
abline(h=mean(dat$npgo.regr, na.rm=T))
abline(v=1989, lty=2)
mtext("l", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,1.5,1.5,1))


# finally, the NPGO-EOF2 plots
lim <- range(npgo.eof2.regr1, npgo.eof2.regr2)
z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern1  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r1  # replace elements NOT corresponding to land with loadings!
z[!land] <- npgo.eof2.regr1  # replace elements NOT corresponding to land with loadings!

z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("m", adj=0.05, line=-1.4, cex=1.1)

z <- rep(NA, ncol(SST))
# z[!land] <- npgo.pattern2  # replace elements NOT corresponding to land with loadings!
# z[!land] <- npgo.eof.r2
z[!land] <- npgo.eof2.regr2
z <- t(matrix(z, length(y.t)))  # Convert vector to matrix and transpose for plotting
image.plot(x.t,y.t,z, col=new.col, zlim=c(lim[1], -lim[1]),
           xlab = "", ylab = "", yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))


contour(x.t,y.t,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', 'Canada', fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires', 'usa',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'USSR',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Japan',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'Mexico',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires', 'China',fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3") 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
mtext("n", adj=0.05, line=-1.4, cex=1.1)

par(mar=c(1.5,3,1.5,1))

plot(dfa.dat$dec.yr, dfa.dat$PDO.NPGO.cor, type="l", xlab="", ylab="PDO-NPGO correlation", xlim=xlim, col=cb[6])
abline(h=0)
abline(v=1989, lty=2)
mtext("o", adj=0.05, line=-1.4, cex=1.1)
dev.off()

###########################
# now plot loadings/trend for dfa

png("basin environment dfa.png", 8, 4, units="in", res=300) 
ggarrange(lp, tp, nrow=1, ncol=2, widths = c(0.5,1))
dev.off()