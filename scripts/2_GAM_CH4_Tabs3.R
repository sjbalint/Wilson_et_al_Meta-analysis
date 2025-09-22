# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson, last update 09.22.25 by Emily Wilson
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Generalized additive model for soil salinity data (GAM), Tab S3 and Fig 2

# -------------------------------------------------------------------------
rm(list=ls()) #clear the environment
#dev.off() #clear plotting window

# install packages --------------------------------------------------------
library(mgcv) #for gam
library(scales) #for transformation function
library(tidyverse) #for data manipulation

# Get the GHG flux data ---------------------------------------------------
#call up the dataset
#meta_data<-readRDS("pw_data.RDS")
meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv")

# Full model --------------------------------------------------------------
#all models weighted by vi (aka n or the number of measurements in each mean, the paper number is a random factor)
gam_model_1a <- gam(yi ~ tide, data = meta_data, weights = vi)
summary(gam_model_1a) #R2=0.31
gam_model_1b <- gam(yi ~ tide + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_1b) #R2=0.31

gam_model_2a <- gam(yi ~ season, data = meta_data, weights = vi)
summary(gam_model_2a) #R2=0.02
gam_model_2b <- gam(yi ~ season + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_2b) #R2=0.08

gam_model_3a <- gam(yi ~ ch4_method_simple, data = meta_data, weights = vi)
summary(gam_model_3a) #R2=0.00
gam_model_3b <- gam(yi ~ ch4_method_simple + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_3b) #R2=0.05

gam_model_4a <- gam(yi ~ salinity, data = meta_data, weights = vi)
summary(gam_model_4a) #R2=0.33
gam_model_4b <- gam(yi ~ salinity + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_4b) #R2=0.33

gam_model_6ab <- gam(yi ~ climate_region, data = meta_data, weights = vi)
summary(gam_model_6ab) #R2=0.441
gam_model_6bb <- gam(yi ~ climate_region + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_6bb) #R2=0.442

#top predictors only
gam_model_5 <- gam(yi ~ s(salinity_conductivity) + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_5) #R2=0.38
gam_model_6 <- gam(yi ~ s(lat) + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_6) #R2=0.52
gam_model_7 <- gam(yi ~ plant_species + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_7) #R2=0.61, 12 plant species, 7 are significant
gam_model_8 <- gam(yi ~ s(lat) + plant_species + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_8) #R2=0.62
gam_model_9 <- gam(yi ~ season + plant_species + s(lat) + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_9) #R2=0.66
gam_model_10 <- gam(yi ~ s(salinity_conductivity) + plant_species +  s(lat) + season + s(p_num, bs="re"), data = meta_data, weights = vi)
summary(gam_model_10) #R2=0.68, 70%

anova(gam_model_9,gam_model_10) #additional predictors are significant
anova(gam_model_8,gam_model_10) #additional predictors are significant

#get all the R2 and dev explained back
summary(gam_model_5)$r.sq 
summary(gam_model_5)$dev
summary(gam_model_6)$r.sq 
summary(gam_model_6)$dev
summary(gam_model_7)$r.sq 
summary(gam_model_7)$dev
summary(gam_model_8)$r.sq 
summary(gam_model_8)$dev 
summary(gam_model_9)$r.sq 
summary(gam_model_9)$dev
summary(gam_model_10)$r.sq 
summary(gam_model_10)$dev

#delta AIC
AIC(gam_model_5)-AIC(gam_model_10)
AIC(gam_model_6)-AIC(gam_model_10)
AIC(gam_model_7)-AIC(gam_model_10)
AIC(gam_model_8)-AIC(gam_model_10)
AIC(gam_model_9)-AIC(gam_model_10)
AIC(gam_model_10)-AIC(gam_model_10)



# Check that the final model meets gam assumptions ------------------------
par(mfrow = c(2, 2))
plot(residuals(gam_model_10)) #normal residuals
qqnorm(residuals(gam_model_10))
hist(residuals(gam_model_10))
dev.off()

par(mfrow = c(2, 2))
gam.check(gam_model_10)
dev.off()

concurvity(gam_model_10)

#convergence in 13 iterations with a very small RMS gradient
#hessian positive definite, no issues with identifiability or optimization.
#model rank: 38/38, efficient use of parameters and no rank deficiency.

#double check most important parameters (cannot weight)
library(gam.hp)
gam_model_10 <- gam(yi ~ s(salinity_conductivity) + plant_species + s(lat) + season + s(p_num, bs="re"), data = meta_data)
summary(gam_model_10) #0.69

# Evaluate relative importance
hp_results <- gam.hp(gam_model_10)

# View the results
print(hp_results) #plants, than lat, than salinity.


#s(salinity_conductivity)  0.0064        0.1217     0.1281     19.37
#plant_species             0.0675        0.2152     0.2827     42.75
#s(lat)                    0.0001        0.2139     0.2140     32.36
#season                    0.0330        0.0033     0.0363      5.49

