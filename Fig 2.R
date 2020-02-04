library(ncdf4)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(oce)
library(R.utils)

# load monthly NCEP/NCAR SLP data!

nc.slp <- nc_open("data/fig2.slp.nc")

# now process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "TIME")  # seconds since 1-1-1970
d <- dates(raw, origin = c(1,1,0001))
m <- months(d)
yr <- years(d)
dec.yr <- as.numeric(as.character(yr)) + (as.numeric(m)-0.5)/12

# and lat/long
x <- ncvar_get(nc.slp, "LON53_101")
y <- ncvar_get(nc.slp, "LAT45_65")

# and extract SLP
SLP <- ncvar_get(nc.slp, "SLP", verbose = F)

# Change data from a 3-D array to a matrix of monthly data
SLP <- aperm(SLP, 3:1)  
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3])) # months in rows, cells in columns! 

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# load pdo
names <- read.table("data/pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("data/pdo", skip=31, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# load npgo
npgo <- read.table("data/npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# limit SLP to NDJ
SLP.NDJ <- SLP[m %in% c("Nov", "Dec", "Jan"),]

yr <- as.numeric(as.character(yr))

# set winter year to correspond to January
win.yr <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr)
# and restrict to Nov-Jan
win.yr <- win.yr[m %in% c("Nov", "Dec", "Jan")]

# get NDJ means for each cell 
rownames(SLP.NDJ) # winter 1949-2019 are complete!

# function to get annual mean
ff <- function(x) tapply(x, win.yr, mean)

# now apply to each cell (column)
SLP.NDJ <- apply(SLP.NDJ, 2, ff)

# limit pdo and npgo to FMA
pdo <- pdo %>%
  filter(month %in% c("FEB", "MAR", "APR"))

PDO.FMA <- tapply(pdo$value, pdo$YEAR, mean) # complete through 2018

npgo <- npgo %>%
  filter(month %in% 2:4)

NPGO.FMA <- tapply(npgo$value, npgo$Year, mean) # complete through 2018

# separate SLP, PDO, and NPGO into era-specific chunks for regression maps

SLP1 <- SLP.NDJ[rownames(SLP.NDJ) %in% 1950:1988,]
SLP2 <- SLP.NDJ[rownames(SLP.NDJ) %in% 1989:2012,]

PDO1 <- PDO.FMA[names(PDO.FMA) %in% 1950:1988]
PDO2 <- PDO.FMA[names(PDO.FMA) %in% 1989:2012]

NPGO1 <- NPGO.FMA[names(NPGO.FMA) %in% 1950:1988]
NPGO2 <- NPGO.FMA[names(NPGO.FMA) %in% 1989:2012]

# calculate separate regressions in each era!
# make objects to catch results
pdo.regr1 <- pdo.regr2 <- npgo.regr1 <- npgo.regr2 <- NA

# now loop through each cell
for(i in 1:ncol(SLP1)){
  #  i <- 1
  mod <- lm(SLP1[,i] ~ PDO1)
  pdo.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP2[,i] ~ PDO2)
  pdo.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP1[,i] ~ NPGO1)
  npgo.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ NPGO2)
  npgo.regr2[i] <- summary(mod)$coef[2,1] 
}


# calculate era differences for each cell
pdo.diff <- pdo.regr2 - pdo.regr1
npgo.diff <- npgo.regr2 - npgo.regr1
diff.lim <- range(pdo.diff, npgo.diff) # limit for plotting

# and a combined plot
png("figs/Fig 2.png", 8,6, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

new.col <- oce.colorsPalette(64)

par(mar=c(1.25,0.5,1,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,3), cex.axis=0.8, cex.lab=0.8, oma=c(0,2,2,0))

xlim <- c(160,250)

lim <- range(pdo.regr1, pdo.regr2, npgo.regr1, npgo.regr2)

# PDO first era
z <- pdo.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico'), 
    fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
 
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=0.5, resolution = 0)
mtext("a) SLP vs. PDO 1950-1988", adj=0)

# PDO second era
z <- pdo.regr2  
z <- t(matrix(z, length(y)))  
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico'), 
    fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=0.5, resolution = 0)
mtext("b) SLP vs. PDO 1989-2012", adj=0)

# PDO diff
z <- pdo.diff 
z <- t(matrix(z, length(y)))
image.plot(x,y,z, col=new.col, zlim=c(-diff.lim[2], diff.lim[2]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico'), 
    fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=0.5, resolution = 0)
mtext("c) SLP vs. PDO (difference)", adj=0)

# npgo first era
z <- npgo.regr1  
z <- t(matrix(z, length(y))) 

image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico'), 
    fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=0.5, resolution = 0)

mtext("d) SLP vs. NPGO 1950-1988", adj=0)

z <- npgo.regr2 
z <- t(matrix(z, length(y)))

image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico'), 
    fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=0.5, resolution = 0)

mtext("e) SLP vs. NPGO 1989-2012", adj=0)

# npgo diff
z <- npgo.diff
z <- t(matrix(z, length(y)))
image.plot(x,y,z, col=new.col, zlim=c(-diff.lim[2], diff.lim[2]), ylim=c(20,68), xlim=xlim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico'), 
    fill=T,xlim=c(130,250), ylim=c(20,70),add=T, lwd=0.5, col="darkgoldenrod3")
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=0.5, resolution = 0)

mtext("f) SLP vs. NPGO (difference)", adj=0)
dev.off()

