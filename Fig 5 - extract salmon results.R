#==================================================================================================
#Project Name: FATE Working Group - Plotting and Compiling Output
#Date: 6.25.19
#
#Purpose: Compile output objects with values for 
#
#
#==================================================================================================
#NOTES:
# 
#
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
library(shinystan)
require(dplyr)

# CONTROL ==========================================================
read <- TRUE #Whether to read in all model files

# Select Species 
species <- c("Sockeye","Pink","Chum")
n.species <- length(species)


# Whether to fit a model with PDO or NPGO 
vars <- c("pdo2a", "npgo2a")
n.vars <- length(vars)

regions <- c("South","GOA","EBS")
n.regions <- length(regions)

# Create Output Objects ==========================================
# By Model
list.rhat <- matrix(nrow=0, ncol=3) # Convergence diag
list.neff <- matrix(nrow=0, ncol=3) # Effective Sample Size

# By Region
list.ratio <- matrix(nrow=0, ncol=4)

# By Stock
list.phi <- matrix(nrow=0, ncol=5) # Autoregressive Coefficient
list.beta <- matrix(nrow=0, ncol=5)
list.beta2 <- matrix(nrow=0, ncol=5)
list.beta_ratio <- matrix(nrow=0, ncol=5)

list.both_beta <- matrix(nrow=0, ncol=6)


# Load A Sample Dataset ==========================================
if(read==TRUE) {

for(s in 1:n.species) {
  s <- 1
  print(paste("s:",s,"of",n.species))
  temp.species <- fit.species <- species[s]
  

  for(v in 1:n.vars) {
    v <- 1
    print(paste("v:",v,"of",n.vars))
    temp.var <- var <- vars[v]
    
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
    
    # Region Lookup
    lookup.stock.region <- data.frame(stocks, stock.regions[region], stringsAsFactors=FALSE)
    
    # Load data ========================
    temp.fit <- readRDS(file=paste0("output/salmon output",temp.species,"-",temp.var,".rds"))
    temp.sum <- summary(temp.fit)
    
    temp.pars <- extract(temp.fit)
    # MODEL: Extract Convergence Diagnostics ==========================
    rhat <- temp.sum$summary[,ncol(temp.sum$summary)]
    neff <- temp.sum$summary[,ncol(temp.sum$summary)-1]
    
    list.rhat <- rbind(list.rhat, data.frame(temp.species, temp.var, rhat))
    list.neff <- rbind(list.neff, data.frame(temp.species, temp.var, neff))
    
    # REGION: Ratio Parameter =============================
    temp.ratio <- data.frame(temp.pars$ratio)
    names(temp.ratio) <- stock.regions
    temp.ratio <- melt(temp.ratio)
    list.ratio <- rbind(list.ratio, data.frame(temp.species, temp.var, temp.ratio))
    
    
    # Extract Autoregression Coeff ===========================================
    temp.phi <- data.frame(temp.pars$phi_trans)
    names(temp.phi) <- stocks
    temp.phi <- melt(temp.phi)
    temp.phi <- temp.phi %>% left_join(lookup.stock.region, by=c("variable"="stocks"))
    list.phi <- rbind(list.phi, data.frame(temp.species, temp.var, temp.phi))
    
    # Beta Coefficient ===========================================
    temp.beta <- data.frame(temp.pars$beta)
    names(temp.beta) <- stocks
    temp.beta <- melt(temp.beta)
    temp.beta <- temp.beta %>% left_join(lookup.stock.region, by=c("variable"="stocks"))
    list.beta <- rbind(list.beta, data.frame(temp.species, temp.var, temp.beta))
    
    temp.beta2 <- data.frame(temp.pars$beta2)
    names(temp.beta2) <- stocks
    temp.beta2 <- melt(temp.beta2)
    temp.beta2 <- temp.beta2 %>% left_join(lookup.stock.region, by=c("variable"="stocks"))
    list.beta2 <- rbind(list.beta2, data.frame(temp.species, temp.var, temp.beta2))
    
    # Both Betas ==========================================
    # both_beta1 <- temp.beta
    # names(both_beta1)[2] <- "Pre"
    # both_beta2 <- temp.beta2
    # names(both_beta2)[2] <- "Post"
    # temp.both_beta <- both_beta1 %>% left_join(both_beta2) %>% gather(key='Period', value="value" -)
    # list.both_beta <- rbind(list.both_beta, temp.both_beta)
    # 
    # Beta Ratio ===========================================
    temp.beta_ratio <- data.frame(temp.pars$beta_ratio)
    names(temp.beta_ratio) <- stocks
    temp.beta_ratio <- melt(temp.beta_ratio)
    temp.beta_ratio <- temp.beta_ratio %>% left_join(lookup.stock.region, by=c("variable"="stocks"))
    list.beta_ratio <- rbind(list.beta_ratio, data.frame(temp.species, temp.var, temp.beta_ratio))
    
    # Traceplots =========================================
    pdf(paste0("figs/salmon diagnostic figs/Traceplot ", temp.species,"_",temp.var,".pdf"), height=10, width=12)
    # plot(rstan::traceplot(temp.fit, pars="mu_beta"))
    # plot(rstan::traceplot(temp.fit, pars="sigma_beta"))
    plot(rstan::traceplot(temp.fit, pars="beta"))
    plot(rstan::traceplot(temp.fit, pars="ratio"))
    plot(rstan::traceplot(temp.fit, pars="ricker_alpha"))
    plot(rstan::traceplot(temp.fit, pars="ricker_beta"))
    plot(rstan::traceplot(temp.fit, pars="sigma_resid"))
    plot(rstan::traceplot(temp.fit, pars="phi"))
    dev.off()
    
  }
}
  
# Save Extracted Objects ==========================================
# MODEL
names(list.neff) <- c("species","var","value")
write.csv(list.neff, "output/salmon outputlist.neff.csv")

names(list.rhat) <- c("species","var","value")
write.csv(list.rhat, "output/salmon outputlist.rhat.csv")

# REGION
names(list.ratio) <- c("species","var","region","value")
write.csv(list.ratio, "output/salmon outputlist.ratio.csv") 

# STOCK
names(list.phi) <- c("species","var","stock","value","region")
write.csv(list.phi, "output/salmon outputlist.phi.csv")

names(list.beta) <- c("species","var","stock","value","region")
write.csv(list.beta, "output/salmon outputlist.beta.csv")

names(list.beta2) <- c("species","var","stock","value","region")
write.csv(list.beta2, "output/salmon outputlist.beta2.csv")

names(list.beta_ratio) <- c("species","var","stock","value","region")
write.csv(list.beta_ratio, "output/salmon outputlist.beta_ratio.csv")

# names(list.mu.ratio) <- c("species","var","region","value")
# names(list.sigma.ratio) <- c("species","var","region","value")
# 
# write.csv(list.mu.ratio, file=file.path(dir.output,"list.mu.ratio.csv")) 
# write.csv(list.sigma.ratio, file=file.path(dir.output,"list.sigma.ratio.csv")) 

}else {
  list.neff <- read.csv("output/salmon outputlist.neff.csv")
  list.rhat <- read.csv("output/salmon outputlist.rhat.csv")
  list.phi <- read.csv("output/salmon outputlist.phi.csv")
  
  list.beta <- read.csv("output/salmon outputlist.beta.csv")
  list.beta2 <- read.csv("output/salmon outputlist.beta2.csv")
  list.ratio <- read.csv("output/salmon outputlist.ratio.csv")
  list.beta_ratio <- read.csv("output/salmon outputlist.beta_ratio.csv")
}
  

# Explore Models with Shiny Stan =======================================
# shiny.fit <- readRDS(file=file.path(dir.output, paste0("Chum-npgo2a.rds")))

# launch_shinystan(shiny.fit) #Required internet
# deploy_shinystan(as.shinystan(shiny.fit))
  
# Plot: rhat ===========================================================
g.rhat <- ggplot(list.rhat, aes(value, fill=var)) +
  theme_linedraw() + 
  geom_density(alpha=0.5) +
  facet_wrap(~species) +
  coord_cartesian(xlim=c(1,1.1))
g.rhat

ggsave(file="figs/salmon diagnostic figs/Rhat.pdf", plot=g.rhat,
       height=4, width=7, units="in")

# Plot: neff ===========================================================
g.neff <- ggplot(list.neff, aes(value, fill=var)) +
  theme_linedraw() + 
  geom_density(alpha=0.5) +
  facet_wrap(~species)
g.neff
ggsave(file="figs/salmon diagnostic figs/Effective Sample Size.pdf", plot=g.neff,
       height=4, width=7, units="in")

# Plot: autocorr ===========================================================
# g.ar <- ggplot(list.phi, aes(value, fill=var)) +
#           theme_linedraw() +
#           geom_density(alpha=0.5) +
#           facet_grid(var~species)
# g.ar

# Plot: ratio ===========================================================

g.ratio <- ggplot(list.ratio, aes(x=region, y=value, fill=var)) +
             scale_fill_viridis_d() +
             theme_linedraw() +
             geom_hline(yintercept=1, lty=1, col='red') +
             geom_violin(alpha=0.5) + 
             facet_wrap(~species) +
             # coord_flip(ylim=c(0,max(g.lims$upper)))
             coord_flip() +
             ylab("Multiplier for Covariate Effect")

g.ratio
g.ratio2 <- g.ratio + facet_grid(species~var)
g.ratio2
ggsave("figs/salmon diagnostic figs/Regional Ratio Param.pdf", plot=g.ratio2,
         height=7, width=8, units="in")



# Beta Ratio ===============================
head(list.beta_ratio)

list.betaRatio <- vector('list', length=n.species)

list.species <- species

for(s in 1:n.species) {
  
  list.betaRatio[[s]] <- list.beta_ratio %>% filter(species==list.species[s]) %>% arrange(region) %>% 
    ggplot(aes(x=stock, y=value, fill=region)) +
    theme_linedraw() +
    scale_fill_viridis_d() +
    geom_hline(yintercept = 0) +
    geom_violin() +
    coord_flip() +
    facet_wrap(~species)
  # ggsave(file=file.path(dir.figs/salmon diagnostic figs, paste0(".pdf")), plot=list.betaRatio[[s]],
  #        height=7, width=8, units="in")
  
} #next s
# Plot combined plots
comb.ratio <- cowplot::plot_grid(plotlist=list.betaRatio, ncol=n.species)
ggsave(file="figs/salmon diagnostic figs/Beta Ratio.pdf", plot=comb.ratio,
       height=7, width=10, units="in")




# Traceplots ===============
rstan::traceplot(temp.fit, pars="beta")










