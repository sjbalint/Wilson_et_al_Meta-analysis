# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson, last update 09.13.25 by Emily Wilson
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Medians and Carbon Offset calculations for Table 3
dev.off()
rm(list=ls()) #clear the environment

library(mgcv) #for gam
library(tidyverse) #for data manipulation
library(scales) #for trans_new

#data_ch4<-readRDS("full_data.RDS") #Full aggregated dataset (n=1031)
data_ch4<-read.csv("raw/ch4_full_dataset.csv")


# Check data distribution -------------------------------------------------
ggplot(data_ch4, aes(x=yi)) +
  geom_histogram(aes(y=after_stat(density)), bins=20, fill="lightblue", color="black") +
  stat_function(fun=dnorm, args=list(mean=mean(data_ch4$yi), sd=sd(data_ch4$yi)), color="red", linewidth=1.2) +
  labs(title="Log CH4") +
  theme_minimal() #looks normal
shapiro.test(data_ch4$yi)  #normal
qqnorm(data_ch4$yi);qqline(data_ch4$yi, col = 2) #normal

#custom asinh transform for ggplot created by Sawyer Balint
asinh_trans <- trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh
)

# Poffenbarger et al. (2011) and VCS predicitions Emmer et al. (2023)---------------------------------------
#For soil >18 GHG=0.011 t CH4 ha-1 yr-1 x GWP
#For soil >20 GHG=0.0056 t CH4 ha-1 yr-1 x GWP
#to get from metric tonnes ha-1 yr-1 to umol m-2 yr-1
#tonnes to grams = 1000000, #ha to m2 = 1/10000, #grams to mols=1/16.04, moles to micromoles = 1000000

#Verra standards
VCS_18<-(0.011*1000000*1000000)/(16.04*10000)
VCS_18 #this is similar to Poffenbarger et al. (2011) at a salinity of 24.5
VCS_20 <-(0.0056*1000000*1000000)/(16.04*10000)
VCS_20 #this is similar to Poffenbarger et al. (2011) at a salinity of 30

#Poffenbarger
#log(ch4) = -0.056xsalinity+1.38, input the salinity into the equation
#10^((-0.056xsalinity)+1.38), untransform
#for grams to mols (1/16.04), for mols to umol (1000000/1)

#predict for salinity categories ---------------------------------------
#Define salinity values, min and max are defined as the highest and lowest salinity within a category
sal_cat_vals <- c(0.5, 2.25, 5,  #oligohaline
                  6, 12, 18,     #mesohaline
                  19, 24.5, 30,  #polyhaline
                  31, 33, 35)    #euhaline

#create conversion for untransformed umol/m2/yr
predict_poff <- function(salinity) {
  (10^((-0.056 * salinity) + 1.38)) * (1 / 16.04) * (1000000 / 1)
}

#apply function to each salinity value
predicted_sal_cat_vals <- data.frame(salinity = sal_cat_vals, ch4_flux = predict_poff(sal_cat_vals))

#check that values look correct
plot(predicted_sal_cat_vals$ch4_flux~predicted_sal_cat_vals$salinity) #shows log-linear relationship from Poffenbarger et al. (2011)!


#Medians -------------------------------------------------------------
full_data <- data_ch4 %>% #remove plant species that only have a few entries
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
seasonal_data <- full_data %>%
  filter(season != "annual") %>% #remove annual values because these are redundant
  group_by(plant_species, salinity, season) %>% #seasonal medians for each plant species and salinity
  summarise(
    n = sum(ch4_n), #how many individual measurements go into the median
    seasonal_median_yi = median(yi, na.rm = TRUE), #get median
    seasonal_q1 = quantile(yi, 0.25),
    seasonal_q3 = quantile(yi, 0.75),
    .groups = "drop"
  )

#create a grid of all the plant and salinity combinations for each season
plants_salinity <- expand_grid(
  plant_species = unique(seasonal_data$plant_species),
  salinity = unique(seasonal_data$salinity),
  season = unique(seasonal_data$season)
)

all_seasonal_data <- plants_salinity %>% #add the seasonal data and remove plants that do not have enough data for specific salinities
  left_join(seasonal_data, by = c("plant_species", "salinity", "season"))
  
write.csv(all_seasonal_data, "output/seasonal_medians.csv")

all_seasonal_data <- all_seasonal_data %>%
  filter(
    !(plant_species == "Cyperus malaccensis" & salinity %in% c("polyhaline", "euhaline")),
    !(plant_species == "Distichlis Spicata" & salinity == "oligohaline"),
    !(plant_species == "Suaeda salsa" & salinity %in% c("polyhaline", "euhaline")),
    !(plant_species == "Scirpus mariqueter" & salinity %in% c("euhaline", "polyhaline")),
    !(plant_species == "Spartina patens" & salinity %in% c("oligohaline", "euhaline")),
    !(plant_species == "Juncus sp." & salinity %in% c("oligohaline", "euhaline")),
    !(plant_species == "Phragmites australis" & salinity == "euhaline"),
    !(plant_species == "Distichlis Spicata")
  ) %>%
  droplevels()


#calculate missing seasonal values using a linear model
lm_model <- lm(seasonal_median_yi ~ salinity + season + plant_species,
               data = all_seasonal_data, na.action = na.exclude)
summary(lm_model) #R2=0.71

q1_model <- lm(seasonal_q1 ~ salinity + season + plant_species,
               data = all_seasonal_data, na.action = na.exclude)
summary(q1_model) #R2=0.57

q3_model <- lm(seasonal_q3 ~ salinity + season + plant_species,
               data = all_seasonal_data, na.action = na.exclude)
summary(q3_model) #R2=0.66

#check which rows are missing the median and Q1 and Q3 and replace with predicted values
all_seasonal_data <- all_seasonal_data %>%
  mutate(predicted_yi = predict(lm_model, newdata = all_seasonal_data)) %>%
  mutate(seasonal_median_yi = ifelse(is.na(seasonal_median_yi), predicted_yi, seasonal_median_yi)) %>%
  select(-predicted_yi)

all_seasonal_data <- all_seasonal_data %>%
  mutate(predicted_q1 = predict(q1_model, newdata = all_seasonal_data)) %>%
  mutate(seasonal_q1 = ifelse(is.na(seasonal_q1), predicted_q1, seasonal_q1)) %>%
  select(-predicted_q1)

all_seasonal_data <- all_seasonal_data %>%
  mutate(predicted_q3 = predict(q3_model, newdata = all_seasonal_data)) %>%
  mutate(seasonal_q3 = ifelse(is.na(seasonal_q3), predicted_q3, seasonal_q3)) %>%
  select(-predicted_q3)

#Calculate annual medians. Multiply by how many hours in a season then add up each season
annual_medians <- all_seasonal_data %>%
  group_by(plant_species, salinity) %>%
  summarise(
    n=sum(n, na.rm = TRUE),
    annual_median_actual = (sum(sinh(seasonal_median_yi)*2190)), # sum of seasonal medians,
    annual_q1_actual = (sum(sinh(seasonal_q1)*2190)), # sum of seasonal medians,
    annual_q3_actual = (sum(sinh(seasonal_q3)*2190)), # sum of seasonal medians
    .groups = "drop"
  )

# Define the salinity categories to add
poff_data <- data.frame(
  plant_species = "Salinity Model",
  salinity = c("oligohaline", "mesohaline", "polyhaline", "euhaline"),
  annual_median_actual = c(1118911.24, 318269.95, 63503.20, 21222.46)
)

annual_medians <- bind_rows(annual_medians, poff_data)

annual_medians<-annual_medians%>% #convert to g/m2/yr
  mutate(annual_median = (annual_median_actual/1000000)*16.04,
         annual_q1 = (annual_q1_actual/1000000)*16.04,
         annual_q3 = (annual_q3_actual/1000000)*16.04,)


# Offset Calculation ------------------------------------------------------
#5.35 Mg C02eq ha year, gas emissions in gas mass units
annual_medians<-annual_medians%>%
  mutate(offset_20 = ((annual_median*96)/ #multiply the the median flux in gas mass units times SGWP and divide by CO2 sequestration
                              (5.35*(1000000/1)*(1/10000))) #convert from MG*ha^-1 to g*m^-2
                            *100) #make it a percent

annual_medians<-annual_medians%>%
  mutate(offset_100 = ((annual_median*45)/ #multiply the the median flux in gas mass units times SGWP and divide by CO2 sequestration
                         (5.35*(1000000/1)*(1/10000))) #convert from MG*ha^-1 to g*m^-2
         *100) #make it a percent

annual_medians <- annual_medians %>%
  mutate(perc_diff = case_when(
    salinity == "oligohaline" ~ ((annual_median-17.94733629) / 17.94733629) * 100,
    salinity == "mesohaline"  ~ ((annual_median-5.10505000) / 5.10505000) * 100,
    salinity == "polyhaline"  ~ ((annual_median-1.018591) / 1.018591) * 100,
    salinity == "euhaline"    ~ ((annual_median-0.34040826) / 0.34040826) * 100,
    TRUE ~ NA_real_
  ))

annual_medians <- annual_medians %>%
  mutate(annual_median = round(annual_median, 2),
         annual_q1 =round(annual_q1, 2),
         annual_q3 =round(annual_q3, 2),
         offset_20 = round(offset_20, 2),
         offset_100 = round(offset_100, 2),
         perc_diff = round(perc_diff, 2))%>%
  select(-annual_median_actual)%>%
  select(-annual_q1_actual)%>%
  select(-annual_q3_actual)
write.csv(annual_medians, "output/annual_medians.csv")

