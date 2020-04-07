#==================================================================================================
#Project Name: FATE Working Group - Hierarchcial Salmon Breakpoint Regression - Free AR parameter (phi)
#Date: 8.23.19
#
#Purpose: Fit hierarchical stock-recruitment model to salmon data that estimates a breakpoint in PDO and NPGO effects
#  a) Fit to each species separately
#  b) Estimate PDO/NPGO effects pre/post breakpoint
#  c) Population level effects assumed to arise from common distribution across locations, but with
#       region-specific mean and standard deviation.
#
#
#==================================================================================================
#NOTES:
# Version 1:
#   Multiplicative effect of era change, SEPARATE PDO/NPGO coefficients (betas) and era change (ratio) variables.
#     Betas and ratios are structured hierarchically, with normal prior the mean and sd of which is region-specific.
# 
#  Version 2:
#    Multiplicative effect of era change, SEPARATE PDO/NPGO coefficients (betas), but era change (ratio) variables
#      are common within regions.
#    Betas (NOT ratios) are structured hierarchically, with normal prior the mean and sd of which is region-specific.
#  
# Version 3:
#    ADDITIVE effect of era change, SEPARATE PDO/NPGO coefficients (betas), but ADDITIVE era change (ratio) variables
#      are common within regions.
#    Betas (NOT ratios) are structured hierarchically, with normal prior the mean and sd of which is region-specific.
# 
# Version 4:
#    ADDITIVE effect of era change, SEPARATE PDO/NPGO coefficients (betas), but ADDITIVE era change (ratio) variables
#      are common within regions.
#    NEITHER betas nor ratios are structured hierarchically. Difference from V3 - Removed hierarchical structure

# Version 5:
#    MULTIPLICATIVE effect of era change, SEPARATE PDO/NPGO coefficients (betas), but MULTIPLICATIVE era change (ratio) variables
#      are common within regions.
#    NEITHER betas nor ratios are structured hierarchically. Difference from V3 - Removed hierarchical structure. REMOVED THE INIT RESIDUAL added autoregressive effect to likelihood. 

# Version 6:
#    MULTIPLICATIVE effect of era change, SEPARATE PDO/NPGO coefficients (betas), but MULTIPLICATIVE era change (ratio) variables
#      are common within regions.
#    NEITHER betas nor ratios are structured hierarchically. Difference from V3 - Removed hierarchical structure. REMOVED THE INIT RESIDUAL added autoregressive effect to likelihood. 
#     Priors: for beta~dnorm(0,1), ratio centered at 1 ratio~dnorm(1,1)


# TIMINGS =====================================
# Version 5:
# [1] "n.iter: 50000"
# [1] "n.thin: 15"
# [1] "Sun Nov  3 11:02:06 2019"
# [1] "Sun Nov  3 13:36:03 2019"

#==================================================================================================
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
library(cowplot)
library(bayesplot)
# library(brms)

dir.create("output/salmon output", showWarnings = F)
dir.create("figs/salmon diagnostic figs", showWarnings = F)

# CONTROL ==========================================================
fit <- TRUE # Do we fit the model, or just load saved .rds outputs

# MCMC Parameters
n.chains <- 3
n.iter <- 5e4
n.thin <- 15
(n.iter/n.thin)*0.5*n.chains

# Select Species
species <- c("Sockeye","Pink","Chum")
fit.species <- species[1]


vars <- c("pdo2a","npgo2a")

var <- vars[1]

# START LOOP =======================================================
# Comment me out if you want to run a single species and variable combination!

# Output Objects
out.waic <- array(dim=c(length(species), length(vars)), dimnames=list(species, vars))
out.se_waic <- array(dim=c(length(species), length(vars)), dimnames=list(species, vars))
out.looic <- array(dim=c(length(species), length(vars)), dimnames=list(species, vars))
out.se_looic <- array(dim=c(length(species), length(vars)), dimnames=list(species, vars))

start <- date()

# fit.species <- species[1]
# var <- vars[1]
for(fit.species in species) {
  for(var in vars) {

   # Define Output File Names ============================================
    file.name <- paste0(fit.species,"-",var)

    # Load/Compile Data ================================================
    # Stock-recruitment data
    dat <- read.csv("data/salmon run dat.csv", header=TRUE, stringsAsFactors=FALSE)
    dat.sr <- dat %>% filter(species==fit.species)
    dat.sr$ln.rps <- log(dat.sr$recruits/dat.sr$spawners)

    # Environmental data
    dat.env <- read.csv("data/winter pdo npgo various smoothings.csv", header=TRUE, stringsAsFactors=FALSE)

    dat.env.fit <- dat.env %>% select(year, var)

    # Bind Environmental data
    dat.input <- dat.sr %>% left_join(dat.env.fit, by=c("entry.yr"="year"))
    # Rename environmental variable
    names(dat.input)[ncol(dat.input)] <- "env.var"

    # Remove any missing environmental observations
    dat.input <- dat.input %>% filter(!is.na(ln.rps) & !is.na(spawners) & !is.nan(ln.rps) & !is.nan(spawners))

    # Metadata
    regions <- unique(dat.input$region)
    n.regions <- length(regions)

    # Create Input Data ================================================
    stocks <- unique(dat.input$stock)
    n.stocks <- length(stocks)

    #Required Attributs of the species
    stock.regions <- unique(dat.input$region)
    n.stock.regions <- length(stock.regions)

    stock.years <- min(dat.input$brood.yr):max(dat.input$brood.yr)
    n.stock.years <- length(stock.years)

    #Create Data Objects
    maxN <- n.stock.years
    S <- n.stocks
    N <- vector(length=n.stocks)
    R <- n.stock.regions #Number of Regions
    region <- vector(length=n.stocks)
    covar <- array(data=0, dim=c(n.stocks,maxN))
    era <- array(data=0, dim=c(n.stocks,maxN))

    #Year Pointer - For Referencing Coefficient Values
    years <- array(data=0, dim=c(n.stocks,maxN))

    #Ricker Parameters
    ln_rps <- array(data=0, dim=c(n.stocks, maxN))
    spawn <- array(data=0, dim=c(n.stocks, maxN))

    p <- 1
    for(p in 1:n.stocks) {

      #Retreive Data =================
      temp.stock <- stocks[p]
      dat.temp <- dat.input %>% filter(stock==temp.stock) %>% arrange(brood.yr)

      #Assign STAN Inputs ============
      N[p] <- nrow(dat.temp) # Number of years in stock-recruitment time series for stock
      region[p] <- which(stock.regions==unique(dat.temp$region)) # Pointer to region for stock

      ln_rps[p, 1:N[p]] <- dat.temp$ln.rps # Response variable for stock
      spawn[p, 1:N[p]] <- dat.temp$spawn # Spawning abundance for density-dependent effect on stock

      years[p,1:N[p]] <- which(stock.years %in% dat.temp$brood.yr ) # Stock-specific brood years

      #Assign Covars and Eras ===============
      covar[p,1:N[p]] <- dat.temp$env.var
      era[p,1:N[p]] <- dat.temp$era

    }#next p


    # Fit Stan Model ===================================================
    #Fit the model
    stan.fit <- NULL
    if(fit==TRUE) {
      stan.fit <- stan(file="models/mod.stan",
                       model_name="mod.stan",
                       data=list("ln_rps"=ln_rps, "spawn"=spawn,
                                 "N"=N, "maxN"=maxN,
                                 "S"=S, "R"=R,
                                 "region"=region,
                                 "covar"=covar, "era"=era
                       ),
                       chains=n.chains, iter=n.iter, thin=n.thin,
                       cores=n.chains, verbose=FALSE,
                       seed=101,
                       control = list(adapt_delta = 0.99))
      #Save Output
      saveRDS(stan.fit, file=paste0("output/salmon output/",file.name,".rds"))
    }else {
      stan.fit <- readRDS(file=paste0("output/salmon output/",file.name,".rds"))
    }
    
    
    # Write a .csv of Model Output =====================================
    write.csv(summary(stan.fit)$summary, file=paste0("output/salmon output/",fit.species, "_",var, "_summary.csv"))
    # Calculate WAIC for Model =========================================


    temp.waic <- waic(extract(stan.fit)$log_lik)
    out.waic[which(species==fit.species), which(vars==var)] <- temp.waic$waic
    out.se_waic[which(species==fit.species), which(vars==var)] <- temp.waic$se_waic

    temp.loo <- loo(extract(stan.fit)$log_lik)
    out.looic[which(species==fit.species), which(vars==var)] <- temp.loo$looic
    out.se_looic[which(species==fit.species), which(vars==var)] <- temp.loo$se_looic

    # END LOOP ========================================================
  } #next var
} #next species

end <- date()

# Save WAIC Looic =================================================
write.csv(out.waic, "output/salmon output/waic.csv")
write.csv(out.se_waic, "output/salmon output/se_waic.csv")
write.csv(out.looic, "output/salmon output/looic.csv")
write.csv(out.se_looic, "output/salmon output/se_looic.csv")


# Timing Diagnostics
print(paste('n.iter:',n.iter))
print(paste('n.thin:',n.thin))
print(start)
print(end)
