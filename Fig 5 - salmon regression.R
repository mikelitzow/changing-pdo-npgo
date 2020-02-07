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

# Define Workflow Paths ============================================
# *Assumes you are working from the Sergent_Streamflow R project
wd <- getwd()
dir.output <- file.path(wd,"output","freeAR_6")
dir.figs <- file.path(wd,"plots","freeAR_6")
dir.data <- file.path(wd,"data")
dir.mods <- file.path(wd, "models")

dir.create(dir.output)
dir.create(dir.figs)

# CONTROL ==========================================================
fit <- TRUE # Do we fit the model, or just load saved .rds outputs

# MCMC Parameters
n.chains <- 3
n.iter <- 5e4#4e4
n.thin <- 15#15
(n.iter/n.thin)*0.5*n.chains

# Select Species
species <- c("Sockeye","Pink","Chum")
fit.species <- species[1]


# Whether to fit a model with PDO or NPGO
# vars <- c("pdo1","pdo2a","pdo2b",
#           "pdo3","npgo","npgo2a",
#           "npgo2b","npgo3")

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
    dat <- read.csv(file.path(dir.data, "salmon run dat.csv"), header=TRUE, stringsAsFactors=FALSE)
    dat.sr <- dat %>% filter(species==fit.species)
    dat.sr$ln.rps <- log(dat.sr$recruits/dat.sr$spawners)

    # Environmental data
    dat.env <- read.csv(file.path(dir.data, "winter pdo npgo various smoothings.csv"), header=TRUE, stringsAsFactors=FALSE)

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


    # Fit brms model ===================================================
    # fit.brm <- brm(ln.rps ~ (1|stock) + spawners:stock + env.var:era + (env.var:era|region/stock),
    #                data=dat.input,
    #                iter=1000, thin=2,
    #                chains=3,
    #                family=gaussian(),
    #                # prior=prior(normal(0,1), class=b))
    #                # prior=c(prior(normal(0,1), class=b),
    #                #         # prior(normal(0,10), class=Intercept),
    #                #         prior(normal(0,10), class=sd),
    #                #         # prior(normal(0,100), class=sds),
    #                #         prior(normal(0,10), class=sds),
    #                #         prior(normal(0,10), class=sigma)),
    #                sample_prior=TRUE,
    #                control = list(adapt_delta = 0.99))
    #
    # summary(fit.brm)

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
    # K <- 2 #Number of covariates PDO, NPGO
    # covars <- array(dim=c(n.stocks,100,K)) #We will start with a temporary length of 100 years then trim down to the max N
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

    #Determine maximum length of covariates =====================
    # maxN <- max(N)
    # temp.regions <- regions[unique(region)]

    # Fit Stan Model ===================================================
    #Fit the model
    stan.fit <- NULL
    if(fit==TRUE) {
      stan.fit <- stan(file=file.path(dir.mods,"hier-Ricker-freeAR_6.stan"),
                       model_name="hier-Ricker-freeAR_6.stan",
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
      saveRDS(stan.fit, file=file.path(dir.output, paste0(file.name,".rds")))
    }else {
      stan.fit <- readRDS(file=file.path(dir.output,paste0(file.name,".rds")))
    }
    
    # temp.plt1 <- stan_trace(stan.fit, pars="ricker_beta")
    # temp.plt2 <- plot(stan.fit, pars='beta')
    
    # Write a .csv of Model Output =====================================
    write.csv(summary(stan.fit)$summary, file=file.path(dir.figs, paste0(fit.species, "_",var, "_summary.csv")))
    # Calculate WAIC for Model =========================================


    temp.waic <- waic(extract(stan.fit)$log_lik)
    out.waic[which(species==fit.species), which(vars==var)] <- temp.waic$waic
    out.se_waic[which(species==fit.species), which(vars==var)] <- temp.waic$se_waic

    temp.loo <- loo(extract(stan.fit)$log_lik)
    out.looic[which(species==fit.species), which(vars==var)] <- temp.loo$looic
    out.se_looic[which(species==fit.species), which(vars==var)] <- temp.loo$se_looic

    # Plot Output ======================================================
    # Someone can continue here.

    # pars <- extract(stan.fit)

    # beta_ratio <- data.frame(pars$mu_ratio)
    # names(mu_ratios) <- stocks
    # exp_mu_ratios <- exp(mu_ratios)
# 
#     list.mu_ratios <- melt(beta_ratio)
# 
#     g <- list.mu_ratios %>% ggplot(aes(value, fill=variable)) +
#       scale_fill_colorblind() +
#       geom_density(alpha=0.5)
#     # g
#     ggsave(file=file.path(dir.figs, paste0(fit.species, "_",var,"mu_ratio hist.png")), plot=g,
#            height=6, width=6, units='in')
# 
#     g2 <- list.mu_ratios %>% ggplot(aes(x=variable, y=value, fill=variable)) +
#       scale_fill_colorblind() +
#       geom_eye(alpha=0.5) +
#       coord_flip() +
#       theme(legend.position = 'none')
#     ggsave(file=file.path(dir.figs,paste0(fit.species, "_",var,"mu_ratio geom_eye.png")), plot=g2,
#              height=6, width=6, units='in')
# 

    # END LOOP ========================================================
  } #next var
} #next species

end <- date()

# Save WAIC Looic =================================================
write.csv(out.waic, file.path(dir.output,"waic.csv"))
write.csv(out.se_waic, file.path(dir.output,"se_waic.csv"))
write.csv(out.looic, file.path(dir.output,"looic.csv"))
write.csv(out.se_looic, file.path(dir.output,"se_looic.csv"))


# Timing Diagnostics
print(paste('n.iter:',n.iter))
print(paste('n.thin:',n.thin))
print(start)
print(end)

























