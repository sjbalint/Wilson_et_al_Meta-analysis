# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson
#last update 09.25.25 by Sawyer Balint
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Generalized additive model for soil salinity data (GAM), Tab S3 and Fig 2

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

#select the best model
#the second model has the same AIC, so when in doubt pick the simplier model
selected.model <- get.models(dredge, 1)[[1]]

summary(selected.model)

# check model assumptions -------------------------------------------------

draw(selected.model)
appraise(selected.model)

# export model ------------------------------------------------------------

saveRDS(selected.model, "Rdata/GAM_SalinityCategory.rds")

