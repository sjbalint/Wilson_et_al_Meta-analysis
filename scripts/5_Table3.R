# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson
#last update 09.29.25 by Sawyer Balint
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Medians and Carbon Offset calculations for Table 3

#rm(list=ls()) #clear the environment

library(mgcv) #for gam
library(tidyverse) #for data manipulation
library(scales) #for trans_new

source("functions/utilities.r")

# import data -------------------------------------------------------------

df <- read.csv("raw/ch4_full_dataset.csv")

# Poffenbarger et al. (2011) and VCS predicitions Emmer et al. (2023)---------------------------------------

#For soil >18 GHG=0.011 t CH4 ha-1 yr-1 x GWP
#For soil >20 GHG=0.0056 t CH4 ha-1 yr-1 x GWP
#to get from metric tonnes ha-1 yr-1 to umol m-2 yr-1
#tonnes to grams = 1000000, #ha to m2 = 1/10000, #grams to mols=1/16.04, moles to micromoles = 1000000

#create conversion for untransformed umol/m2/yr
predict_poff <- function(salinity) {
  log.g.CH4 <- (-0.056 * salinity) + 1.38 #equation from Poffenbarger (2011)
  g.CH4 <- 10^(log.g.CH4) #back-transfoprm
  umol.CH4 <- g.CH4 * (1000000/16.04) #convert from g to umol
  return(umol.CH4)
}

#predict for salinity categories ---------------------------------------

#Define salinity values, min and max are defined as the highest and lowest salinity within a category
salinity.df <- data.frame(salinity=c("oligohaline", "mesohaline", "polyhaline", "euhaline"),
                          min=c(0.5, 6, 19, 31),
                          mean=c(2.25, 12, 24.5, 33),
                          max=c(5,18,30,35)) %>%
  pivot_longer(-salinity, names_to="stat", values_to="salinity.ppt") %>%
  mutate(salinity=factor(salinity, levels=unique(salinity)),
         CH4_poff = predict_poff(salinity.ppt))

ggplot(salinity.df, aes(salinity.ppt, CH4_poff, fill=salinity, shape=salinity))+
  theme_classic()+
  geom_point()+
  scale_shape_manual(values=c(21:26))+
  scale_fill_viridis_d()+
  labs(x="Salinity (ppt)",
       y=bquote(Predicted~CH[4]~Flux~(mu*mol~m^-2~yr^-1)))+
  theme(legend.position="inside",
        legend.position.inside=c(0.8, 0.8),
        legend.title=element_blank())

#Medians -------------------------------------------------------------

subset.df <- df %>% #remove plant species that only have a few entries
  filter(
    !plant_species %in% c(
      "Puccinellia maritima", "Plantago maritima", "Spartina foliosa",
      "Elytrigia repens", "Spartina angelica", "Carex sp.",
      "Schoenoplectus triqueter", "Aeluropus littoralis",
      "Sporobolus virginicus", "Schoenoplectus americanus",
      "Tamarix chinensis", "Cladium jamaicense", "Salicornia sp.", "Suaeda sp."
    )
  )

#calculate seasonal medians and Q1, Q2
seasonal.df <- subset.df %>%
  filter(season != "annual") %>% #remove annual values because these are redundant
  group_by(plant_species, salinity, season) %>% #seasonal medians for each plant species and salinity
  summarise(
    n = sum(ch4_n), #how many individual measurements go into the median
    yi_median = median(yi, na.rm = TRUE), #get median
    yi_q1 = quantile(yi, 0.25),
    yi_q3 = quantile(yi, 0.75),
  ) %>%
  ungroup()

#create a grid of all the plant and salinity combinations for each season
grid.df <- expand_grid(
  plant_species = unique(seasonal.df$plant_species),
  salinity = unique(seasonal.df$salinity),
  season = unique(seasonal.df$season)
)

seasonal.df <- left_join(grid.df, seasonal.df)
  
write.csv(seasonal.df, "output/seasonal_medians.csv")

#remove missing and limited species/salinity combinations
count.df <- seasonal.df %>%
  drop_na(n) %>%
  group_by(plant_species, salinity) %>%
  summarize(n.seasons = length(unique(season)),
            n.total = sum(n, na.rm=TRUE))  %>%
  ungroup() %>%
  filter(n.seasons>1,
         n.total>15)

seasonal.df <- left_join(count.df, seasonal.df)

# model missing seasons ---------------------------------------------------

#calculate missing seasonal values using a linear model
lm_model <- lm(yi_median ~ salinity + season + plant_species,
               data = seasonal.df, na.action = na.exclude)

summary(lm_model) #R2=0.71

q1_model <- lm(yi_q1 ~ salinity + season + plant_species,
               data = seasonal.df, na.action = na.exclude)

summary(q1_model) #R2=0.57

q3_model <- lm(yi_q3 ~ salinity + season + plant_species,
               data = seasonal.df, na.action = na.exclude)

summary(q3_model) #R2=0.67

# fill in missing data with the linear models -----------------------------

#check which rows are missing the median and Q1 and Q3 and replace with predicted values
seasonal.df <- seasonal.df %>%
  mutate(yi_median_pred = predict(lm_model, newdata = .),
         yi_q1_pred = predict(q1_model, newdata = .),
         yi_q3_pred = predict(q3_model, newdata = .),
         yi_median = ifelse(is.na(yi_median), yi_median_pred, yi_median),
         yi_q1 = ifelse(is.na(yi_q1), yi_q1_pred, yi_q1),
         yi_q3 = ifelse(is.na(yi_q3), yi_q3_pred, yi_q3)) %>%
  select(-c(yi_median_pred, yi_q1_pred, yi_q3_pred))

#Calculate annual medians. Multiply by how many hours in a season then add up each season
median.df <- seasonal.df %>%
  group_by(plant_species, salinity) %>%
  summarise(
    n=sum(n, na.rm = TRUE),
    CH4_umolm2yr_median = (sum(sinh(yi_median)*2190)), # sum of seasonal medians, #hours per season
    CH4_umolm2yr_q1 = (sum(sinh(yi_q1)*2190)), # sum of seasonal medians,
    CH4_umolm2yr_q3 = (sum(sinh(yi_q3)*2190)), # sum of seasonal medians
  ) %>%
  ungroup()

# Define the salinity categories to add
poff.df <- salinity.df %>%
  filter(stat=="mean") %>%
  select(salinity, CH4_poff)

annual_medians <- poff.df %>%
  mutate(plant_species = "Salinity Model") %>%
  rename(CH4_umolm2yr_median = CH4_poff) %>%
  bind_rows(median.df, .) %>% #convert to g/m2/yr
  mutate(CH4_gm2yr_median = (CH4_umolm2yr_median/1000000)*16.04, #molar mass of CH4
         CH4_gm2yr_q1 = (CH4_umolm2yr_q1/1000000)*16.04,
         CH4_gm2yr_q3 = (CH4_umolm2yr_q3/1000000)*16.04)

# Offset Calculation ------------------------------------------------------

#5.35 Mg C02eq ha year, gas emissions in gas mass units
MgCO2.ha.year = 5.35
gCO2.ha.year = MgCO2.ha.year*1000000 #Mg to g
gCO2.m2.year = gCO2.ha.year/10000 #hectare to m2

#multiply the the median flux in gas mass units times SGWP and divide by CO2 sequestration
offset.df <- left_join(annual_medians, poff.df) %>%
  mutate(CH4_poff = (CH4_poff/1000000)*16.04,
         offset_20 = 100*(CH4_gm2yr_median*96)/gCO2.m2.year,
         offset_100 = 100*(CH4_gm2yr_median*45)/gCO2.m2.year,
         perc_diff = 100 * (CH4_gm2yr_median - CH4_poff) / CH4_poff
  )

# export data -------------------------------------------------------------

output.df <- offset.df %>%
  select(plant_species, salinity, n, CH4_gm2yr_median, CH4_gm2yr_q1, CH4_gm2yr_q3, offset_20, offset_100, perc_diff) %>%
  mutate_if(is.numeric, round, digits=2)

write.csv(output.df, "output/annual_medians.csv")

