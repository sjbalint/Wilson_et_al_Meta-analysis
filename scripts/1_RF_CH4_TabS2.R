# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson, last update 09.22.25 by Emily Wilson
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis
#Random Forest Models for Table S2 and Figure 2

#dev.off() #clear plotting window
rm(list=ls()) #clear the environment

# install packages --------------------------------------------------------
library(tidyverse) #for data manipulation
library(randomForestSRC) #parallel random forest
library(scales) #for transformation function

# Get the GHG flux data ---------------------------------------------------
# With PW conductivity ----------------------------------------------------
#meta_data <- readRDS("pw_data.RDS") #get data
meta_data<-read.csv("raw/ch4_soilsalinity_dataset.csv")

hist(meta_data$yi) #ch4 data is called yi, it was transformed using asinh(yi)
qqnorm(meta_data$yi) #asinh preserves negatives and 0 and makes the data normal

# Prepare final meta_data -------------------------------------------------------
meta_data <- meta_data %>% #select only what we need
  select(yi, lat, plant_species, season, salinity, salinity_conductivity, ch4_method_simple, tide)

meta_data<-data.frame(meta_data)
meta_data<-meta_data%>% #categorical variables need to be factors
  mutate(
    plant_species=as.factor(plant_species),
    salinity=as.factor(salinity),
    tide=as.factor(tide),
    ch4_method_simple=as.factor(ch4_method_simple),
    season=as.factor(season),
    lat=abs(as.numeric(lat)))
str(meta_data)
#tune(yi~., data=meta_data, trace=TRUE)
set.seed(39)
model <- rfsrc(yi~.,
               ntree=5000, 
               data=meta_data, 
               nodesize = 12, #lowest error, highest r2
               mtry=5, #suggest by tune
               case.wt = meta_data$vi,
               forest=TRUE,
               importance = TRUE)
#looks good
plot(model)
print(model)


# Make a Variable Importance Plot -----------------------------------------

var_imp <- data.frame(
  Variable = names(model$importance),
  Importance = model$importance)

var_imp <- var_imp %>%
  mutate(Variable = case_when(
    Variable == "tide" ~ "Tidal Range",
    Variable ==  "climate_region" ~ "Climate Region",
    Variable ==  "season" ~ "Season",
    Variable ==  "salinity" ~ "Salinity Category",
    Variable == "ch4_method_simple" ~ "Sampling Method",
    Variable == "plant_species" ~ "Plant Species",
    Variable == "lat" ~  "Absolute Latitude",
    Variable == "salinity_conductivity"~"Soil salinity",
    TRUE ~   Variable
  ))

#arrange from highest to lowest importance
var_imp <- var_imp %>%
  arrange(desc(Importance))

p2<-ggplot(var_imp, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "darkgreen", color="black") +
  coord_flip() +
  theme_classic()+
  labs(
    x = "Variable",
    y = "Importance")+ 
  theme(
    axis.text.x = element_text(size = 10, colour = "black", angle = 20, hjust = 0.9, family = "Helvetica"), 
    axis.text.y = element_text(size = 12, colour = "black", family = "Helvetica"),
    axis.title.x = element_text(size = 11, color = "black", family = "Helvetica"),
    axis.title.y = element_text(size = 12, color = "black", family = "Helvetica"))
p2
ggsave("figures/CH4_forest.png",width=6,height=5)

