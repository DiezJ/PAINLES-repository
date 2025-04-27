
# Load stuff ----
library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggplot2)
library(ggthemes)
library(modelr)
library(ggExtra)
library(rstan)

library(lme4)
library(lmerTest)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
theme_set(theme_minimal())


# dat.new <- readRDS("Data/dat.new.rds")
# head(dat.new,2); dim(dat.new)
# summary(dat.new$PCA1)
# sapply(dat.new, function(x) sum(is.na(x)))
# 

# dat.leafN_OLD <- readRDS("Data/dat.leafN.rds")
# dat.LDMC_OLD <- readRDS("Data/dat.LDMC.rds")
# dat.SLA_OLD <- readRDS("Data/dat.SLA.rds")

# head(dat.leafN_OLD,1); dim(dat.leafN_OLD)
# head(dat.LDMC_OLD,1); dim(dat.LDMC_OLD)
# head(dat.SLA_OLD,1); dim(dat.SLA_OLD)

# Using "new" larger datasets
dat.leafN <- readRDS("Data/dat.leafN.rds")
dat.LDMC <- readRDS("Data/dat.LDMC.rds")
dat.SLA <- readRDS("Data/dat.SLA.rds")

head(dat.leafN,1); dim(dat.leafN)
head(dat.LDMC,1); dim(dat.LDMC)
head(dat.SLA,1); dim(dat.SLA)


use.regs <- c("Northwest Forests","Deserts","Great Plains","Eastern Forests","Northern Forests") # Dana's order
regs<-use.regs
# did this in processing
# dat.new <- filter(dat.new, ecoregion %in% use.regs)
# dat.leafN <- filter(dat.leafN, ecoregion %in% use.regs)
# head(dat.leafN,2); dim(dat.leafN)
# dat.LDMC <- filter(dat.LDMC, ecoregion %in% use.regs)
# dat.SLA <- filter(dat.SLA, ecoregion %in% use.regs)

dim(dat.leafN)                                                                             
table(dat.leafN$GrowthForm)

dat.leafN <- droplevels(dat.leafN)
table(dat.leafN$GrowthForm, dat.leafN$ExoticStatus,dat.leafN$ecoregion)


## Scaling ----
ln.scaled <- dat.leafN$ln.RelCover_scaled
ln <- dat.leafN$ln.RelCover
RA <- exp(dat.leafN$ln.RelCover)
# plot(ln, ln.scaled)
# plot(reg, ln)
# ln.scaled = (ln - mean(ln)) / sd(ln)
# to go from ln.scaled units back to Rel.abund
# unscaled: ln.scaled * sd(ln) + mean(ln)
ln.calc <- ln.scaled * sd(ln) + mean(ln)
# plot(ln.calc, ln); abline(0,1,col='red')
# to go from ln units back to Rel.abund
# Rel.Abund: exp(ln)
RA_calc <- exp(ln.calc)
RA_calc2 <- exp(ln.scaled * sd(ln) + mean(ln))
# plot(RA_calc, RA); abline(0,1,col='red')
# plot(RA_calc2, RA_calc); abline(0,1,col='red')


table(dat.leafN$GrowthForm)

# export data for Dana
# write.table(dat.LDMC, "dat.LDMC.csv", sep=',', row.names=FALSE, eol = "\r")
# write.table(dat.SLA, "dat.SLA.csv", sep=',', row.names=FALSE, eol = "\r")
# write.table(dat.leafN, "dat.leafN.csv", sep=',', row.names=FALSE, eol = "\r")

## Models: ----

# Leaf N ----

## m3 ----

  m3.leafN.brms <- bf(ln.RelCover_scaled ~  GrowthForm + ExoticStatus*ln.leafN_scaled*CWMnonfocal.LeafN_scaled + (1|SpCode) + (ln.leafN_scaled|Plot))
  
  for(r in 1:length(regs)){
    # r <- 5
    reg <- regs[r]  
    dat_sub <- filter(dat.leafN, ecoregion == reg)
    # dat <- filter(dat.leafN, ecoregion == "Great Plains")
    head(dat_sub,n=2); dim(dat_sub)
    table(dat_sub$GrowthForm)
    
    fitm3 <- brm(
      m3.leafN.brms, data = dat_sub, 
      prior(normal(0,1), class="b"),
      prior(exponential(1), class="sd"),
      prior(exponential(1), class="sigma"),
      family = gaussian(),
      control = list(adapt_delta = 0.9, max_treedepth = 15),
      cores = 4,
      chains=4,
      ,init="0"
      ,refresh=200
      ,iter=4000
    )
    # saveRDS(fitm3, file=paste("Saved Models/",reg,"_leafN_fitm3_RandomEffect.rds", sep=''))
    
    # fitm3
    
    saveRDS(fitm3, file=paste("Saved Models/",reg,"_leafN_fitm3.rds", sep=''))
    # saveRDS(fitm3, file="Saved Models/Great Plains_leafN_fitm3.rds")
    
    
    posterior_cp <- as.array(fitm3)
    pars1 <- dimnames(posterior_cp)[[3]]
    pars1 <- pars1[grep('b_', pars1)]
    np_cp <- nuts_params(fitm3)
    pdf(file=paste("Saved Models/",reg,"_leafN_m3.traceplots.pdf",sep=''), width = 12, height = 10)
    print( mcmc_trace(posterior_cp, pars = pars1, np = np_cp) + 
             xlab("Post-warmup iteration")
    )
    dev.off()
    
  } # end region loop
  

# LDMC & SLA, the same ----

