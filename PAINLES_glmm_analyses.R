
#GLMM analyses of trait differences between native and non-native species

# Packages and options
library(glmmTMB)
library(DHARMa)
library(emmeans)

options(contrasts = c("contr.sum","contr.poly"))


##Load data for regional-scale trait analyses
painles_regional<-read.table("PAINLES_dat8_region.csv",sep=",",head=TRUE,stringsAsFactors = FALSE)
attach(painles_regional)
names(painles_regional)
str(painles_regional)

## Leaf N regional analysis
leaf_N_regional <- glmmTMB(lnLeafN ~ NativeStatus*ecoregion + GrowthForm,
                           na.action=na.omit,
                           data = painles_regional)

summary(leaf_N_regional)
simulationOutput <- simulateResiduals(fittedModel = leaf_N_regional, plot = T)
emmeans(leaf_N_regional, poly ~ NativeStatus | ecoregion)

## SLA regional analysis
SLA_regional <- glmmTMB(lnSLA ~ NativeStatus*ecoregion + GrowthForm,
                        na.action=na.omit,
                        data = painles_regional)

summary(SLA_regional)
simulationOutput <- simulateResiduals(fittedModel = SLA_regional, plot = T)
emmeans(SLA_regional, poly ~ NativeStatus | ecoregion)

## LDMC regional analysis
LDMC_regional <- glmmTMB(lnLDMC ~ NativeStatus*ecoregion + GrowthForm,
                         na.action=na.omit,
                         data = painles_regional)

summary(LDMC_regional)
simulationOutput <- simulateResiduals(fittedModel = LDMC_regional, plot = T)
emmeans(LDMC_regional, poly ~ NativeStatus | ecoregion)


###Plot scale analyses of trait differences between native and non-native species

##Load data for N analysis at the plot scale
painles_plot_N <-read.table("2origins_leafN_v3.csv",sep=",",head=TRUE,stringsAsFactors = FALSE)
attach(painles_plot_N)
names(painles_plot_N)
str(painles_plot_N)

##Plot scale analysis Leaf N
N_plot <- glmmTMB(lnleafN~ NativeStatus*ecoregion + GrowthForm + (1|Plot),
                  na.action=na.omit,
                  data = painles_plot_N)
summary(N_plot)
simulationOutput <- simulateResiduals(fittedModel = N_plot, plot = T)
emmeans(N_plot, poly ~ NativeStatus | ecoregion)

##Load data for SLA analysis at the plot scale
painles_plot_SLA <-read.table("2origins_SLA_v3.csv",sep=",",head=TRUE,stringsAsFactors = FALSE)
attach(painles_plot_SLA)
names(painles_plot_SLA)
str(painles_plot_SLA)

##Plot scale analysis Leaf SLA
SLA_plot <- glmmTMB(lnSLA~ NativeStatus*ecoregion + GrowthForm + (1|Plot),
                    na.action=na.omit,
                    data = painles_plot_SLA)
summary(SLA_plot)
simulationOutput <- simulateResiduals(fittedModel = SLA_plot, plot = T)
emmeans(SLA_plot, poly ~ NativeStatus | ecoregion)

##Load data for LDMC analys1s at the plot scale
painles_plot_LDMC <-read.table("2origins_LDMC_v3.csv",sep=",",head=TRUE,stringsAsFactors = FALSE)
attach(painles_plot_LDMC)
names(painles_plot_LDMC)
str(painles_plot_LDMC)

##Plot scale analysis Leaf LDMC
LDMC_plot <- glmmTMB(lnLDMC~ NativeStatus*ecoregion + GrowthForm + (1|Plot),
                     na.action=na.omit,
                     data = painles_plot_LDMC)
summary(LDMC_plot)
simulationOutput <- simulateResiduals(fittedModel = LDMC_plot, plot = T)
emmeans(LDMC_plot, poly ~ NativeStatus | ecoregion)


##PAINLES GLMM analysis of cover differences between native and non-native species

##Load data for cover analysis
Mean_cov <-read.table("painles_cover_2024_1Yr_V2.csv",sep=",",head=TRUE,stringsAsFactors = FALSE)
attach(Mean_cov)
names(Mean_cov)
str(Mean_cov)

## Mean cover analysis
Mean_cover <- glmmTMB(LnCov ~ NativeStatus*ecoregion + GrowthForm + (1|Plot) + (1|SpCode),
                    na.action=na.omit,
                    data = Mean_cov)

summary(Mean_cover)
simulationOutput <- simulateResiduals(fittedModel = Mean_cover, plot = T)
emmeans(Mean_cover, poly ~ NativeStatus | ecoregion)

