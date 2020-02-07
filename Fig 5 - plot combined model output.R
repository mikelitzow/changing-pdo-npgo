library(tidyverse)
library(ggpubr)

# load salmon results
salmon.ratio <- read.csv("output/freeAR_6/list.ratio.csv", row.names = 1) 

# rename covariates!
salmon.ratio$var <- ifelse(salmon.ratio$var=="pdo2a", "PDO - salmon", "NPGO - salmon")
salmon.ratio$var <- reorder(salmon.ratio$var, desc(salmon.ratio$var))

# rename region as system and reassign factor levels
names(salmon.ratio)[3] <- "system"
salmon.ratio$system <- ifelse(salmon.ratio$system=="EBS", "Bering Sea",
                              ifelse(salmon.ratio$system=="GOA", "Gulf of Alaska", "Northern California Current"))

# # add dummy lines to make shared legend show up correctly?
# dummy.lines <- data.frame(species=NA, var=NA, 
#                           system=c("Central California Current", "Southern California Current"),
#                           value=NA)
# 
# salmon.ratio <- rbind(salmon.ratio, dummy.lines)

# order the systems north-south
salmon.ratio$order <- ifelse(salmon.ratio$system=="Bering Sea", 1,
                         ifelse(salmon.ratio$system=="Gulf of Alaska", 2,
                                ifelse(salmon.ratio$system=="Northern California Current", 3,
                                       ifelse(salmon.ratio$system=="Central California Current", 4, 5))))

salmon.ratio$system <- reorder(salmon.ratio$system, salmon.ratio$order)

# and order the spp
salmon.ratio$spp.order <- ifelse(salmon.ratio$species=="Pink", 1, 
                                 ifelse(salmon.ratio$species=="Sockeye", 2, 3))

salmon.ratio$species <- reorder(salmon.ratio$species, salmon.ratio$spp.order)


# colorblind...
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plotting functions
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.90 <- function(x) { return(quantile(x, probs=c(0.05,0.95))) }

salmon.plt <- ggplot(salmon.ratio, aes(x=reorder(system, desc(system)), y=value, fill=system)) +
  theme_bw() +
  scale_fill_manual(values=cb[c(6,3,4)],
                    labels=c("Bering Sea", "Gulf of Alaska",
                             "Northern Cal. Curr.")) +
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  # stat_summary(fun.y="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun.y="q.90", colour="black", geom="line", lwd=0.5) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_grid(species~var) +
  ylab("Era 2 slope : Era 1 slope") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_line(size=0),
        legend.position = 'none') +
  geom_hline(aes(yintercept=1), size=0.5) +
  coord_flip(ylim=c(-2,2))

salmon.plt


# load environment results
pdo.env.data <- read.csv("output/pdo_environment_model_data.csv", row.names = 1)
npgo.env.data <- read.csv("output/npgo_environment_model_data.csv", row.names = 1)

npgo.env.data$var  <- "NPGO - environment"
pdo.env.data$var <- "PDO - environment"

# combine and plot pdo/npgo and environment data
env.data <- rbind(pdo.env.data, npgo.env.data)

env.data$var.order <- ifelse(env.data$var=="PDO - environment", 1, 2)
env.data$var <- reorder(env.data$var, env.data$var.order)

# order the systems north-south
env.data$order <- ifelse(env.data$system=="Bering Sea", 1,
                         ifelse(env.data$system=="Gulf of Alaska", 2,
                                ifelse(env.data$system=="Northern California Current", 3,
                                       ifelse(env.data$system=="Central California Current", 4, 5))))

env.data$system <- reorder(env.data$system, env.data$order)

env.plt <- ggplot(env.data, aes(x=reorder(system, desc(system)), y=ratio/100, fill=system)) +
  theme_bw() +
  scale_fill_manual(values=cb[c(6,3,4,2,8)],
                    labels=c("Bering Sea", "Gulf of Alaska",
                             "Northern Cal. Curr.", "Central Cal. Curr.", "Southern Cal. Curr.")) +
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.90", colour="black", geom="line", lwd=0.5) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=2) +
  theme(axis.text.y = element_blank(), axis.title = element_blank(), axis.ticks.y = element_line(size=0),
        legend.position = 'none') +
  geom_hline(aes(yintercept=1), size=0.5) +
  coord_flip(ylim=c(-3,3))

env.plt

# and now the biology results

pdo.biol.data <- read.csv("output/pdo_biology_model_data.csv", row.names=1)
npgo.biol.data <- read.csv("output/npgo_biology_model_data.csv", row.names=1)

npgo.biol.data$var  <- "NPGO - other biology"
pdo.biol.data$var <- "PDO - other biology"

# combine and plot pdo/npgo and biolironment data
biol.data <- rbind(pdo.biol.data, npgo.biol.data)

biol.data$var.order <- ifelse(biol.data$var=="PDO - other biology", 1, 2)
biol.data$var <- reorder(biol.data$var, biol.data$var.order)

# order the systems north-south
biol.data$order <- ifelse(biol.data$system=="Bering Sea", 1,
                         ifelse(biol.data$system=="Gulf of Alaska", 2,
                                ifelse(biol.data$system=="Northern California Current", 3,
                                       ifelse(biol.data$system=="Central California Current", 4, 5))))

biol.data$system <- reorder(biol.data$system, biol.data$order)

biol.plt <- ggplot(biol.data, aes(x=reorder(system, desc(system)), y=ratio/100, fill=system)) +
  theme_bw() +
  scale_fill_manual(values=cb[c(6,3,4,2,8)],
                    labels=c(" Bering Sea", " Gulf of Alaska",
                             " Northern Cal. Curr.", " Central Cal. Curr.", " Southern Cal. Curr.")) +
  
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.90", colour="black", geom="line", lwd=0.5) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=2) +
  ylab("Era 2 slope : Era 1 slope") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_line(size=0),
        legend.position = 'right', legend.title = element_blank()) +
  geom_hline(aes(yintercept=1), size=0.5) +
  coord_flip(ylim=c(-3,3))

biol.plt

# save legend!
plot.leg <- get_legend(biol.plt)
as_ggplot(plot.leg)
ggsave("plots/legend.plot.png")

# and redefine without a legend
biol.plt <- ggplot(biol.data, aes(x=reorder(system, desc(system)), y=ratio/100, fill=system)) +
  theme_bw() +
  scale_fill_manual(values=cb[c(6,3,4,2,8)],
                    labels=c("Bering Sea", "Gulf of Alaska",
                             "Northern Cal. Curr.", "Central Cal. Curr.", "Southern Cal. Curr.")) +
  
  geom_violin(alpha = 0.75, lwd=0.1, scale='width') +
  stat_summary(fun.y="q.90", colour="black", geom="line", lwd=0.5) +
  stat_summary(fun.y="q.50", colour="black", geom="line", lwd=1) +
  stat_summary(fun.y="median", colour="black", size=2, geom="point", pch=21) +
  facet_wrap(~var, ncol=2) +
  ylab("Era 2 slope : Era 1 slope") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_line(size=0),
        legend.position = 'none') +
  geom_hline(aes(yintercept=1), size=0.5) +
  coord_flip(ylim=c(-3,3))
# now combine all three
# make a blank plot as a placeholder

null.plot <- ggplot() + theme_void()

png("plots/combined environment salmon other biology bayesian results plot.png",
    8, 8, units="in", res=300)
ggarrange(null.plot, env.plt, salmon.plt, biol.plt,
          ncol=2, nrow=2)
dev.off()

# and now (finally) make a table of the distributions

biol.table <- biol.data %>%
  group_by(system, var) %>%
  summarise(percntl.10=q.90(ratio/100)[1],
            percntl.25=q.50(ratio/100)[1],
            median=median(ratio/100),
            percntl.75=q.50(ratio/100)[2],
            percntl.90=q.90(ratio/100)[2])

biol.table$group <- "Other biology"

biol.table$var <- ifelse(biol.table$var == "PDO - other biology", "PDO", "NPGO")

# environment results
env.table <- env.data %>%
  group_by(system, var) %>%
  summarise(percntl.10=q.90(ratio/100)[1],
            percntl.25=q.50(ratio/100)[1],
            median=median(ratio/100),
            percntl.75=q.50(ratio/100)[2],
            percntl.90=q.90(ratio/100)[2])

env.table$group <- "Environment"

env.table$var <- ifelse(env.table$var == "PDO - environment", "PDO", "NPGO")

# salmon results
salm.table <- salmon.ratio %>%
  group_by(system, var, species) %>%
  summarise(percntl.10=q.90(value)[1],
            percntl.25=q.50(value)[1],
            median=median(value),
            percntl.75=q.50(value)[2],
            percntl.90=q.90(value)[2])

salm.table$group <- ifelse(salm.table$species == "Pink", "Pink salmon",
                           ifelse(salm.table$species == "Sockeye", "Sockeye salmon", "Chum salmon"))
salm.table <- salm.table %>%
  select(-species)

salm.table$var <- ifelse(salm.table$var == "PDO - salmon", "PDO", "NPGO")

SI.table <- rbind(env.table, salm.table, biol.table)

write.csv(SI.table, "output/SI posterior distributions table.csv", row.names = F)
