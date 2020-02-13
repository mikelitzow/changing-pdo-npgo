library(tidyverse)
library(zoo)
library(nlme)
library(MARSS)
library(MuMIn)
library(ggpubr)
library(reshape2)
library(glmmTMB)
library(rstanarm)
library(lme4)
library(rstan)
library(ggthemes)
library(tidybayes)
library(bayesplot)

# download PDO / NPGO and process
# download.file("http://jisao.washington.edu/pdo/PDO.latest", "data/pdo")
names <- read.table("data/pdo", skip=30, nrows=1, as.is = T)
pdo <- read.table("data/pdo", skip=32, nrows=119, fill=T, col.names = names)
pdo$YEAR <- 1900:(1899+nrow(pdo)) # drop asterisks!
pdo <- pdo %>%
  gather(month, value, -YEAR) %>%
  arrange(YEAR)

# download.file("http://www.oces.us/npgo/npgo.php", "~npgo")
npgo <- read.table("data/npgo", skip=10, nrows=828, fill=T, col.names = c("Year", "month", "value"))

# calculate NDJFM means for each index
pdo$win.yr <- ifelse(pdo$month %in% c("NOV", "DEC"), pdo$YEAR+1, pdo$YEAR)

# limit to winter months only
pdo <- pdo %>%
  filter(month %in% c("NOV", "DEC", "JAN", "FEB", "MAR"))

win.pdo <- tapply(pdo$value, pdo$win.yr, mean)

npgo$win.yr <- ifelse(npgo$month %in% 11:12, npgo$Year+1, npgo$Year)

# limit to winter months only
npgo <- npgo %>%
  filter(month %in% c(11,12,1:3))

win.npgo <- tapply(npgo$value, npgo$win.yr, mean)

# and smoothed (2yr) values of each
win.npgo <- rollapply(win.npgo, 2, mean, align="right", fill=NA)
names(win.npgo) <- 1950:2019
win.pdo <- rollapply(win.pdo, 2, mean, align="right", fill=NA)
names(win.pdo) <- 1900:2019

# load five "non-salmon" data sets: EBS groundfish recruitment, GOA crustaceans/fish,Farallon seabirds, CalCOFI ichthyo,
# CCC seabirds
dat <- read.csv("data/farallon.sbrd.biol.csv", row.names = 1)
# add era term
dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m1 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m1$system <- "Central California Current"

####
# CalCOFI
dat <- read.csv("data/calcofi.biol.csv", row.names=1)

# examine distributions
look <- dat %>%
  gather(key, value, -year)

ggplot(look, aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales="free")

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m2 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m2$system <- "Southern California Current"

#########
# Gulf of Alaska
dat <- read.csv("data/goa.biol.csv")
colnames(dat)[1] <- "year"

# examine distributions
look <- dat %>%
  gather(key, value, -year)

ggplot(look, aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales="free")

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
# standardize all the time series by variable -- so slopes are on same scale
m3 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m3$system <- "Gulf of Alaska"

#######
# Eastern Bering Sea
dat <- read.csv("data/ebs.biol.data.csv")
dat[,2:5] <- sqrt(dat[,2:5])

# examine distributions
look <- dat %>%
  gather(key, value, -year)

ggplot(look, aes(value)) +
  geom_histogram() +
  facet_wrap(~key, scales="free")

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

# reshape with year, era, and pdo and npgo as the grouping variables
melted <- melt(dat, id.vars = c("year","pdo","era","npgo"))
melted$variable_era = paste0(melted$era,melted$variable)
melted$value <- as.numeric(melted$value)

# standardize all the time series by variable -- so slopes are on same scale
m4 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m4$system <- "Bering Sea"

#######
# Northern California Current
dat <- read.csv("data/ncc.biol.dat.csv")

# examine distributions (only one non-salmon time series!)
ggplot(dat, aes(value)) +
  geom_histogram() 

dat$era <- as.factor(ifelse(dat$year <= 1988, 1, 2))

# and pdo/npgo
dat$pdo <- win.pdo[match(dat$year, names(win.pdo))]
dat$npgo <- win.npgo[match(dat$year, names(win.npgo))]

melted <- dat
melted$variable_era = paste0(melted$era,melted$variable)
melted$value <- as.numeric(melted$value)

# standardize
m5 = dplyr::group_by(melted, variable) %>%
  mutate(scale_x = scale(value))
m5$system <- "Northern California Current"

#######
# combine

melted <- rbind(m1, m2, m3, m4, m5)
melted$year <- as.numeric(melted$year)
melted$variable <- as.factor(melted$variable)
melted$variable_era <- as.factor(melted$variable_era)


# save a version to add to SI
write.csv(melted, "data/regional.non-salmon.biology.data.for.SI.csv", row.names = F)



# make an object to capture model output
model.data = data.frame()

# set up to loop through each region/ecosystem separately
levels.syst <- as.factor(unique(melted$system))

for(s in levels.syst) {

  temp <- melted %>%
    filter(system==s)
  temp <- na.omit(temp)

  # create data for stan
  temp$variable = as.character(temp$variable)
  temp$variable = as.factor(temp$variable)

  stan_data = list(era = as.numeric(temp$era),
                   y = temp$pdo,
                   variable = as.numeric(temp$variable),
                   n = nrow(temp),
                   n_levels = max(as.numeric(temp$variable)),
                   x = temp$scale_x)

  mod = stan(file="models/mod.stan", data=stan_data, chains=3, warmup=4000, iter=6000,thin=2,
    pars = c("beta","mu_beta","ratio","mu_ratio","sigma_beta","sigma_ratio",
      "exp_mu_ratio","exp_ratio","pred"),
    control=list(adapt_delta=0.99, max_treedepth=20))

  pars = rstan::extract(mod,permuted=TRUE)

  model.data = rbind(model.data,
    data.frame(system=s,
      ratio=100*(pars$mu_ratio)))

  temp$pred = apply(pars$pred,2,mean)

  # create an array for caterpillar plots by variable
  draws = rstan::extract(mod,permuted=FALSE)
  par_names = dimnames(draws)$parameters
  par_names[grep("^mu_ratio", par_names)] = "global mean"
  par_names[grep("^ratio", par_names)] = levels(temp$variable)
  dimnames(draws)$parameters = par_names
  idx = which(par_names %in% c("global mean",levels(temp$variable)))

}

# order the systems north-south
model.data$order <- ifelse(model.data$system=="Bering Sea", 1,
                           ifelse(model.data$system=="Gulf of Alaska", 2,
                                  ifelse(model.data$system=="Northern California Current", 3,
                                         ifelse(model.data$system=="Central California Current", 4, 5))))

model.data$system <- reorder(model.data$system, model.data$order)

pdo.biol.data <- model.data

# save for future reference
write.csv(pdo.biol.data, "output/pdo_biology_model_data.csv")

#################
## and the same thing for npgo
model.data <- data.frame()
for(s in levels.syst) {

  temp <- melted %>%
    filter(system==s)
  temp <- na.omit(temp)
  # create data for stan
  temp$variable = as.character(temp$variable)
  temp$variable = as.factor(temp$variable)

  stan_data = list(era = as.numeric(temp$era),
                   y = temp$npgo,
                   variable = as.numeric(temp$variable),
                   n = nrow(temp),
                   n_levels = max(as.numeric(temp$variable)),
                   x = temp$scale_x)

  mod = stan(file="models/mod.stan", data=stan_data, chains=3, warmup=4000, iter=6000,thin=2,
    pars = c("beta","mu_beta","ratio","mu_ratio","sigma_beta","sigma_ratio",
      "exp_mu_ratio","exp_ratio","pred"),
    control=list(adapt_delta=0.99, max_treedepth=20))

  pars = rstan::extract(mod,permuted=TRUE)

  model.data = rbind(model.data,
    data.frame(system=s, ratio=100*(pars$mu_ratio)))

  # temp$pred = apply(pars$pred,2,mean)

  # create an array for caterpillar plots by variable
  draws = rstan::extract(mod,permuted=FALSE)
  par_names = dimnames(draws)$parameters
  par_names[grep("^mu_ratio", par_names)] = "global mean"
  par_names[grep("^ratio", par_names)] = levels(temp$variable)
  dimnames(draws)$parameters = par_names
  idx = which(par_names %in% c("global mean",levels(temp$variable)))

}


# order the systems north-south
model.data$order <- ifelse(model.data$system=="Bering Sea", 1,
                           ifelse(model.data$system=="Gulf of Alaska", 2,
                                  ifelse(model.data$system=="Northern California Current", 3,
                                         ifelse(model.data$system=="Central California Current", 4, 5))))

model.data$system <- reorder(model.data$system, model.data$order)

npgo.biol.data <- model.data

# save for future reference
write.csv(npgo.biol.data, "output/npgo_biology_model_data.csv")
