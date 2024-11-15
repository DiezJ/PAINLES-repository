
# Packages ----
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

# obtain and read in raw data ----
# dat.leafN <- readRDS("Data/dat.leafN.rds")
# dat.LDMC <- readRDS("Data/dat.LDMC.rds")
# dat.SLA <- readRDS("Data/dat.SLA.rds")

use.regs <- c("Northwest Forests","Deserts","Great Plains","Eastern Forests","Northern Forests") 
regs<-use.regs

dat.leafN <- filter(dat.leafN, ecoregion %in% use.regs)
dat.LDMC <- filter(dat.LDMC, ecoregion %in% use.regs)
dat.SLA <- filter(dat.SLA, ecoregion %in% use.regs)
dat.leafN <- droplevels(dat.leafN)

# Models: ----
# Leaf N ----

## m3: model with trait data ----

  m3.leafN.brms <- bf(ln.RelCover_scaled ~  GrowthForm + ExoticStatus*ln.leafN_scaled*CWMnonfocal.LeafN_scaled + (1|SpCode) + (ln.leafN_scaled|Plot))

# Running model separately for each region (for computational reasons; could in theory be all run together)
for(r in 1:length(regs)){
  reg <- regs[r]  
  dat_sub <- filter(dat.leafN, ecoregion == reg)

  fitm3 <- brm(
    m3.leafN.brms, data = dat_sub, 
    prior = c( prior(normal(0,1), class="b")),
    control = list(adapt_delta = 0.9, max_treedepth = 15),
    cores = 3,
    chains=3,
    ,init="0"
    ,refresh=100
    ,iter=5000
  )
  
  # saveRDS(fitm3, file=paste("Saved Models/",reg,"_leafN_fitm3.rds", sep=''))
  
# get and save traceplots and model diagnostics
posterior_cp <- as.array(fitm3)
pars1 <- dimnames(posterior_cp)[[3]]
pars1 <- pars1[grep('b_', pars1)]
np_cp <- nuts_params(fitm3)
# pdf(file=paste("Saved Models/",reg,"_leafN_m3.traceplots.pdf",sep=''), width = 12, height = 10)
# print( mcmc_trace(posterior_cp, pars = pars1, np = np_cp) + 
#          xlab("Post-warmup iteration")
# )
# dev.off()
}

## m0 - traitless model ----

# m0.brms <- bf(ln.RelCover_scaled ~  GrowthForm + ExoticStatus + (1 | Plot) + (1|SpCode)) # main one have been using
 m0.brms <- bf(ln.RelCover_scaled ~  GrowthForm + ExoticStatus + (1 | Plot) + (1|SpCode) - 1)
# m0.brms <- bf(ln.RelCover_scaled ~  (1 | GrowthForm) + ExoticStatus + (1 | Plot) + (1|SpCode)) # Growth form as RE

for(r in 1:length(regs)){
  
  dat.region <- droplevels(filter(dat.leafN, ecoregion == regs[r])) # using leaf N dataset
  
  fitm0 <- brm(
    m0.brms, data = dat.region, 
    prior = c( prior(normal(0,1), class="b")),
    control = list(adapt_delta = 0.9, max_treedepth = 15),
    cores = 3,
    chains=3,
    # thin=2
    ,init="0"
    ,refresh=200
    # backend = "cmdstanr", threads = threading(3),
    # set_cmdstan_path(path = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/cmdstanr"),
    ,iter=5000
  )
  
  saveRDS(fitm0, file=paste("Saved Models/",regs[r],"_fitm0_leafN_RandomEffect_laptop.rds",sep=''))

    saveRDS(fitm0, file=paste("Saved Models/",regs[r],"_fitm0_leafN_Minus1.rds",sep=''))


    saveRDS(fitm0, file=paste("Saved Models/",regs[r],"_fitm0_leafN.rds",sep=''))
  
  
  posterior_cp <- as.array(fitm0)
  pars1 <- dimnames(posterior_cp)[[3]][grep( "b_", dimnames(posterior_cp)[[3]])]
  # np_cp <- nuts_params(fitm0)
  pdf(file=paste("Saved Models/",regs[r],"_fitm0_leafN.rds.diagnost.pdf",sep=''), width = 12, height = 10)
  print( mcmc_trace(posterior_cp, pars = pars1) + # , np = np_cp
           xlab("Post-warmup iteration")
  )
  # plot.new()
  # grid.table(fitm0.print)
  dev.off()
}


# Note: Identical models used for other leaf traits (LDMC, SLA)


