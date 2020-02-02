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
weights <-  sqrt(cos(lat*pi/180))
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

# and calculate standard deviation over 11-year rolling windows
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

ggsave("AL temporal variance.png", width=3, height=2, units="in")
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

# and plot
png("figs/Fig 1.png", 4, 7, units="in", res=300) 
ggarrange(a.plot, b.plot, c.plot, labels = c("a)", "b)", "c)"),  nrow=3, align="v")
dev.off()
