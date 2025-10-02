# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson, last update 09.26.25 by Emily Wilson and Sawyer Balint
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Random Forest Models for Table S2 and Figure 2

#dev.off() #clear plotting window
rm(list=ls()) #clear the environment

# install packages --------------------------------------------------------
library(tidyverse) #for data manipulation
library(randomForestSRC) #parallel random forest
library(scales) #for transformation function

# import data -------------------------------------------------------------

df <- read.csv("raw/ch4_soilsalinity_dataset.csv")

#check data distribution
hist(df$yi) #ch4 data is called yi, it was transformed using asinh(yi)
qqnorm(df$yi) #asinh preserves negatives and 0 and makes the data normal

# Prepare final meta_data -------------------------------------------------------

input.df <- df %>%
  select(yi, lat, plant_species, season, salinity, salinity_conductivity, ch4_method_simple, tide) %>%
  #categorical variables need to be factors
  mutate(
    plant_species=as.factor(plant_species),
    salinity=as.factor(salinity),
    tide=as.factor(tide),
    ch4_method_simple=as.factor(ch4_method_simple),
    season=as.factor(season),
    lat=abs(as.numeric(lat)))

str(input.df)

#tune(yi~., data=meta_data, trace=TRUE)

set.seed(39)

model <- rfsrc(yi~.,
               ntree=5000, 
               data=input.df, 
               nodesize = 12, #lowest error, highest r2
               mtry=5, #suggest by tune
               case.wt = input.df$vi,
               forest=TRUE,
               importance = TRUE)

#check model validity
plot(model)
print(model)

# export model ------------------------------------------------------------

saveRDS(model, "Rdata/RandomForest.rds")

