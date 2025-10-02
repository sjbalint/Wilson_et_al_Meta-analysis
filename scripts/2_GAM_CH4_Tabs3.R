# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson
#last update 09.26.25 by Emily Wilson and Sawyer Balint
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Generalized additive model for soil salinity data (GAM), Tab S3, Tab S4, Tab S5, and Fig 2

# -------------------------------------------------------------------------
rm(list=ls()) #clear the environment

# install packages --------------------------------------------------------
library(mgcv) #for gam
library(scales) #for transformation function
library(tidyverse) #for data manipulation
library(gam.hp) 
library(gratia) #for gam diagnostics
library(MuMIn) #for model selection

# import data -------------------------------------------------------------

meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv") %>%
  mutate(plant_species = factor(plant_species),
         salinity = factor(salinity),
         season = factor(season))

# model of soil conductivity ----------------------------------------------

global.model <- gam(yi ~ 
                      s(salinity_conductivity) +
                      plant_species + 
                      s(lat) + 
                      season + 
                      s(p_num, bs="re"), 
                    na.action=na.fail, #needed for dredge
                    data = meta_data, 
                    method="REML",
                    weights = vi)

dredge <- dredge(global.model)

#best model via AIC does not include random effects
dredge

all_models <- get.models(dredge, subset = TRUE)  # gets all models

model_stats <- lapply(all_models, function(m) {
  s <- summary(m)
  data.frame(
    R2 = s$r.sq,
    dev_explained = s$dev.expl
  )
})

# combine into a single data.frame
model_stats_df <- do.call(rbind, model_stats)

dd_df <- as.data.frame(dredge)
dd_stats_salinity <- cbind(dd_df, model_stats_df)

#select the best model
selected.model <- get.models(dredge, 1)[[1]]

summary(selected.model)

# check model assumptions -------------------------------------------------

draw(selected.model)
appraise(selected.model)

# model of soil conductivity ----------------------------------------------
#salinity by plant species

global.model <- gam(yi ~ 
                    s(salinity_conductivity, by=plant_species) +
                    plant_species + 
                    s(lat) + 
                    season + 
                    s(p_num, bs="re"), 
                  na.action=na.fail, #needed for dredge
                  data = meta_data, 
                  method="REML",
                  weights = vi)

dredge <- dredge(global.model)

#best model via AIC does not include random effects
dredge

all_models <- get.models(dredge, subset = TRUE)  # gets all models

model_stats <- lapply(all_models, function(m) {
  s <- summary(m)
  data.frame(
    R2 = s$r.sq,
    dev_explained = s$dev.expl
  )
})

# combine into a single data.frame
model_stats_df <- do.call(rbind, model_stats)

dd_df <- as.data.frame(dredge)
dd_stats_salinitybyplants <- cbind(dd_df, model_stats_df)

#select the best model including random effects
selected.model <- get.models(dredge, 1)[[1]]

summary(selected.model)

# check model assumptions -------------------------------------------------

draw(selected.model)
appraise(selected.model)

# export model ------------------------------------------------------------

saveRDS(selected.model, "Rdata/GAM_SoilSalinity.rds")

# model of salinity category ----------------------------------------------

global.model <- gam(yi ~ 
                      salinity+
                      plant_species+
                      salinity:plant_species +
                      s(lat) + 
                      season + 
                      s(p_num, bs="re"), 
                    na.action=na.fail, #needed for dredge
                    data = meta_data, 
                    method="REML",
                    weights = vi)

# model selection ---------------------------------------------------------

dredge <- dredge(global.model)

#best model via AIC does not include random effects
dredge

all_models <- get.models(dredge, subset = TRUE)  # gets all models

model_stats <- lapply(all_models, function(m) {
  s <- summary(m)
  data.frame(
    R2 = s$r.sq,
    dev_explained = s$dev.expl
  )
})

# combine into a single data.frame
model_stats_df <- do.call(rbind, model_stats)

dd_df <- as.data.frame(dredge)
dd_stats_salinity_category <- cbind(dd_df, model_stats_df)

#select the best model
selected.model_1 <- get.models(dredge, 1)[[1]]
selected.model_2 <- get.models(dredge, 3)[[1]]
anova(selected.model_1,selected.model_2)
summary(selected.model_2)

# check model assumptions -------------------------------------------------

draw(selected.model_2)
appraise(selected.model_2)

# export model ------------------------------------------------------------

saveRDS(selected.model_2, "Rdata/GAM_SalinityCategory.rds")

