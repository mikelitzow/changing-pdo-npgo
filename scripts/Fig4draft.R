library(tidyverse)
library(mgcv)

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# first, pdo/npgo node plots
dat <- read.csv("SKT_Classes.csv", na.strings = "FALSE")

head(dat)

dat <- dat %>%
  arrange(Year)

dat <- na.omit(dat)

dat$order <- ifelse(dat$index=="PDO", 1, 2)
dat$index <- reorder(dat$index, dat$order)

dat$node <- as.factor(dat$node)

dat$era <- ifelse(dat$Year<=1988, 1, 2)
dat$era.label <- as.factor(ifelse(dat$era==1, "Before 1988/89", "After 1988/89"))
dat$era.label <- reorder(dat$era.label, dat$era)

png("pgo npgo node TS.png",6,4, units="in", res=300)
ggplot(dat, aes(node, value)) +
  theme_linedraw() + 
  geom_text(aes(label=Year, color=era.label), size=3) + 
  facet_wrap(~index) + 
  scale_color_manual(values=cb[c(6,2)]) +
  theme(legend.title = element_blank(), legend.position = "top")

dev.off()

# SST time series plot
# load SST
# load ERSSTv5 data
# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2019-04-01)][(0.0):1:(0.0)][(30):1:(66)][(150):1:(250)]", "~updated.sst")
# uncomment that line to download the data

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

# calculate time series

f <- function(x) tapply(x, m[yr %in% 1951:1980], mean)  # function to compute monthly means for a single time series

mu <- apply(SST[yr %in% 1951:1980,], 2, f)	# Compute monthly means for each time series (location)

mu <- mu[rep(1:12, floor(length(d.t)/12)),] 

xtra <- 12*((length(d.t)/12)-floor(length(d.t)/12))

mu <- rbind(mu, mu[1:xtra,])

weights <-  sqrt(cos(lat.t*pi/180))
ff <- function(x) weighted.mean(x, w=weights, na.rm=T)

SST.anom <- SST-mu
SST.anomTS <- apply(SST.anom,1,ff)   # Compute monthly anomalies!

dec.yr.t <- as.numeric(as.character(yr)) + (as.numeric(m)-0.5)/12

# sst time series plot
use.dat <- data.frame(yr=dec.yr.t, sst=SST.anomTS)
use.dat$era <- as.factor(ifelse(use.dat$yr<1989,1,2))

pred.sst <- predict(gam(sst ~ s(yr), data=use.dat))

png("ne pacific sst plot.png", 5,4, units="in", res=300)
par(mar=c(4,4,1.5,1),las=1)
plot(dec.yr.t, SST.anomTS, type="l", xlab="", ylab="ºC wrt 1951-1980", col=cb[6])
lines(dec.yr.t, pred.sst, col="red")
abline(h=0)
abline(v=1989, lty=2)
dev.off()

####
# now pdo-npgo correlations

# load TS
# load pdo
names <- read.table("~pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("~pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load npgo
# download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("~npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))
pdoTS <-pdo[601:1425,]
npgoTS <- npgo[1:825,]
rownames(pdoTS) <- rownames(npgoTS) <- 1:nrow(pdoTS)
npgo$dec.yr <- npgo$Year + (npgo$month-0.5)/12

# create object to catch results and loop through
pdo.npgo.plot <- data.frame(dec.yr=NA, cor=NA)

for(i in 127:(nrow(pdoTS)-126)){
  temp <- data.frame(dec.yr=NA, cor=NA)
  
  pdo.npgo.plot$cor[(i-126)] <- cor(pdoTS$value[(i-126):(i+126)], npgoTS$value[(i-126):(i+126)])
  pdo.npgo.plot$dec.yr[(i-126)] <- npgo$dec.yr[i]

  pdo.npgo.plot<- rbind(pdo.npgo.plot, temp)  
}

png("pdo npgo corr plot.png", 5,4, units="in", res=300)
par(mar=c(4,4,1.5,1),las=1)
plot(pdo.npgo.plot$dec.yr, pdo.npgo.plot$cor, type="n", xlab="", ylab="Correlation")
abline(h=0)
abline(v=1989, lty=2)
lines(pdo.npgo.plot$dec.yr, pdo.npgo.plot$cor,  col=cb[6])
dev.off()