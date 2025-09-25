# -------------------------------------------------------------------------
#Created 1.13.25 by Emily Wilson
#last update 09.25.25 by Sawyer Balint
#Title: Plant Species Drive Global Coastal Wetland Methane Fluxes Meta-analysis

#Supplemental Figure
rm(list=ls()) #clear the environment

# install packages --------------------------------------------------------
library(scales) #for stats
library(tidyverse) #for data manipulation
library(patchwork) #to put plots together
library(purrr)
library(broom)

source("functions/utilities.r")

# import data -------------------------------------------------------------

df <- read.csv("raw/ch4_soilsalinity_dataset.csv")

gam.model <- readRDS("Rdata/GAM_SoilSalinity.rds")

# Plot salinity and yi -------------------------------------------------------

meta_data <- df %>%
  filter(!is.na(salinity_conductivity), #no winter bc Poffenbarger et al. is growing season
         season!="winter") %>% #clean up the data
  filter(!is.na(salinity_conductivity),
         salinity_conductivity < 60) %>%  #only include fluxes with soil salinity
  mutate(lat=abs(as.numeric(lat)), #make lat absolute for analysis
         plant_species=as.factor(plant_species),
         season=as.factor(season),
         salinity_conductivity=round(as.numeric(salinity_conductivity), 2),
         ch4_method_simple=as.factor(ch4_method_simple),
         tide=as.factor(tide))%>%
  filter(plant_species!="Juncae sp")


# create sequence of x values
line_data <- data.frame(x = seq(0,
                                max(meta_data$salinity_conductivity, na.rm = TRUE),
                                length.out = 100)) %>%
  mutate(y = (((-0.056 * x) + 1.38))/0.34,
         y = 10^y,
         y = y*(1/24)*(1/16.04)*(1000/1))

# figure S4 ---------------------------------------------------------------

#get r2

model1<-lm(yi~salinity_conductivity, data=meta_data)

r2 <- summary(model1)$r.squared

r2_label <- bquote("This Study's Relationship between Salinity and "*CH[4])

r2_label <- expression("_ This Study's Relationship between Salinity and CH"[4] * " flux, R"^2 * " = 0.10")

# Get predictions with prediction interval
preds <- predict(model1, newdata = meta_data, interval = "prediction", se.fit = TRUE)
meta_data <- meta_data %>%
  mutate(
    fit = exp(preds$fit[, "fit"]),
    fit.low = exp(preds$fit[, "lwr"]),
    fit.high = exp(preds$fit[, "upr"]))

p1<-ggplot(meta_data, aes(x = salinity_conductivity)) + 
  geom_point(aes(y = sinh(yi), shape=plant_species, fill=plant_species), size=2, alpha = 0.6) +
  geom_ribbon(aes(ymin = fit.low, ymax = fit.high), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit), color = "black", linewidth = 1) +
  geom_vline(xintercept=18, linetype="dashed", color="darkblue", linewidth=1)+
  geom_hline(yintercept=0, linetype="solid", color="black")+
  geom_line(data = line_data, aes(x = x, y =y), color = "darkred", linetype="dotdash", linewidth=1) +
  scale_y_continuous(transform = asinh_trans,
                     breaks=c(-1000,-100,-10,0,10,100,1000,10000))+
  theme_classic()+
  scale_shape_manual(values = c("Cyperus malaccensis"=21, "Spartina alterniflora"=22,
                                "Distichlis Spicata" = 8, "Spartina patens" = 24, 
                                "Suaeda salsa" = 25, "Phragmites australis" = 23,
                                "Cladium jamaicense" = 13, "Juncus sp." = 21,
                                "Juncae sp" = 21,
                                "Sesuvium portulacastrum" = 18, "Plantago maritima" = 10,
                                "Scirpus mariqueter" = 11, "Tamarix chinensis"=14,
                                "Salicornia sp."=12, "Suaeda sp." = 3))+
  scale_fill_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                               "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                               "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                               "Cladium jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#30BCAD",
                               "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black")) +
  theme(legend.position="right")+
  theme(legend.title=element_blank())+
  annotate("text", x = 31, y = 2200,  # Adjust x and y based on your plot range
           label = expression(paste(". _ . _ Salinity-Based Proxy Equation Line")),
           size = 2.3, color = "darkred", family = "Arial")+
  annotate("text", x = 31, y = 3500,  # Adjust x and y based on your plot range
           label = expression(paste("---- 18 ppt Salinity Threshold")),
           size = 2.3, color = "darkblue", family = "Arial")+
  annotate("text", x = 30, y = 6000, label = r2_label, 
           , size = 2.3, family = "Arial" ) +
  theme(
    legend.text =element_text(size = 7, colour = "black", family = "Arial"),
    axis.text.x = element_text(size = 7, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
  labs(x = "Salinity (ppt)", fill = "Plant Species", shape="Plant Species", y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), family="Arial")

p1




meta_data<-df %>%
  filter(!is.na(salinity_conductivity), #no winter bc Poffenbarger et al. is growing season
         season!="winter") %>% #clean up the data
  filter(!is.na(salinity_conductivity),
         salinity_conductivity < 60) %>%  #only include fluxes with soil salinity
  mutate(lat=abs(as.numeric(lat)), #make lat absolute for analysis
         plant_species=as.factor(plant_species),
         season=as.factor(season),
         salinity_conductivity=round(as.numeric(salinity_conductivity), 2),
         ch4_method_simple=as.factor(ch4_method_simple),
         tide=as.factor(tide))%>%
  filter(plant_species!="Juncae sp")

summary(lm(yi~salinity_conductivity, data=meta_data))

meta_data <- meta_data %>%
  mutate(repeat_flag = ifelse("repeat" == "Y", TRUE, FALSE)) %>%
  mutate(group_key = ifelse(!repeat_flag,
                            paste(p_num, lat, ch4_chamber, estuary, marsh, plant_species, salinity, specific_site, section, sep = "_"),
                            NA)) %>%
  # Create numeric unique IDs for each group_key
  mutate(unique_id = ifelse(!repeat_flag,
                            as.integer(factor(group_key)),
                            NA)) %>%
  relocate(unique_id, .after = p_num)

meta_data <- meta_data %>%
  group_by(unique_id) %>%
  summarise(
    # Calculate mean of yi
    yi = mean(yi, na.rm = TRUE),
    
    # Calculate sum of vi (n = number of rows per group)
    vi = as.numeric(sum(vi)), #used to weight gams
    
    ch4_n=as.integer(sum(ch4_n)),
    
    salinity_conductivity = round(as.numeric(mean(salinity_conductivity), 2)),
    
    #Aggregate the columns for data analysis
    p_num = as.factor(paste(unique(p_num), collapse = ", ")),
    plant_species = as.factor(paste(unique(plant_species), collapse = ", ")),
    salinity = as.factor(paste(unique(salinity), collapse = ", ")),
    marsh = as.factor(paste(unique(marsh), collapse = ", ")),
    season = as.factor(paste(unique(season), collapse = ", ")),
    tide = as.factor(paste(unique(tide), collapse = ", ")),
    lat = paste(unique(lat), collapse = ", "),  
    tide = as.factor(paste(unique(tide), collapse = ", ")),
    ch4_method_simple = as.factor(paste(unique(ch4_method_simple), collapse = ", ")),
    climate_region = as.factor(paste(unique(climate_region), collapse = ", ")),
    month = as.factor(paste(unique(month), collapse = ", ")),
    #sampling_year = as.numeric(paste(unique(sampling_year), collapse = ", "))
    
    # Add other columns as needed here...
  )

# Fit model
model1 <- lm(yi ~ salinity_conductivity, data = meta_data)

# Get R²
r2 <- summary(model1)$r.squared
r2
r2_label <- expression("_ This Study's Relationship between Salinity and CH"[4] * " flux, R"^2 * " = 0.06")


# Get predictions with prediction interval
preds <- predict(model1, newdata = meta_data, interval = "prediction", se.fit = TRUE)
meta_data <- meta_data %>%
  mutate(
    fit = exp(preds$fit[, "fit"]),
    fit.low = exp(preds$fit[, "lwr"]),
    fit.high = exp(preds$fit[, "upr"]))

p2<-ggplot(meta_data, aes(x = salinity_conductivity)) + 
  geom_point(aes(y = sinh(yi), shape=plant_species, fill=plant_species), size=2, alpha = 0.6) +
  geom_ribbon(aes(ymin = fit.low, ymax = fit.high), alpha = 0.1, fill = "black") +
  geom_line(aes(y = fit), color = "black", linewidth = 1) +
  geom_vline(xintercept=18, linetype="dashed", color="darkblue", linewidth=1)+
  geom_hline(yintercept=0, linetype="solid", color="black")+
  geom_line(data = line_data, aes(x = x, y =y), color = "darkred", linetype="dotdash", linewidth=1) +
  scale_y_continuous(transform = asinh_trans,
                     breaks=c(-1000,-100,-10,0,10,100,1000,10000))+
  theme_classic()+
  scale_shape_manual(values = c("Cyperus malaccensis"=21, "Spartina alterniflora"=22,
                                "Distichlis Spicata" = 8, "Spartina patens" = 24, 
                                "Suaeda salsa" = 25, "Phragmites australis" = 23,
                                "Cladium jamaicense" = 13, "Juncus sp." = 21,
                                "Juncae sp" = 21,
                                "Sesuvium portulacastrum" = 18, "Plantago maritima" = 10,
                                "Scirpus mariqueter" = 11, "Tamarix chinensis"=14,
                                "Salicornia sp."=12, "Suaeda sp." = 3))+
  scale_fill_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                               "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                               "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                               "Cladium jamaicense" =  "#564e95", "Juncus sp." = "#86a291", "Juncae sp" = "lightblue", 
                               "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#30BCAD",
                               "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black")) +
  theme(legend.position="none")+
  annotate("text", x = 31, y = 2200,  # Adjust x and y based on your plot range
           label = expression(paste(". _ . _ Salinity-Based Proxy Equation Line")),
           size = 2.3, color = "darkred", family = "Arial")+
  annotate("text", x = 31, y = 3500,  # Adjust x and y based on your plot range
           label = expression(paste("---- 18 ppt Salinity Threshold")),
           size = 2.3, color = "darkblue", family = "Arial")+
  annotate("text", x = 30, y = 6000, label = r2_label, 
           , size = 2.3, family = "Arial" ) +
  theme(
    axis.text.x = element_text(size = 7, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 7, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 7, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 7, color = "black", family = "Arial")) +
  labs(x = "Salinity (ppt)", fill = "Plant Species", shape="Plant Species", y = expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")), family="Arial")
p2
p2 <- p2 + theme(legend.position = "none")
p2 <- p2 + guides(fill = "none", color = "none", shape = "none")


combined_plot <- (p1/p2) +
  plot_layout(nrow = 2,guides = "collect") +
  plot_annotation(
    tag_levels = 'a', 
    tag_suffix = "."  # use lower-case letters
  ) &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    plot.tag = element_text(family = "Arial Black", size = 20)
  )

combined_plot

ggsave("figures/FigureS4.png", combined_plot, width = 170, height = 260, dpi=300, units = "mm")

# figure S5 ---------------------------------------------------------------

meta_data <- df %>%
  group_by(plant_species) %>%
  filter(n() > 15) %>%
  ungroup()

newdata.df <- expand_grid(plant_species = unique(meta_data$plant_species),
                          salinity_conductivity = seq(
                            min(meta_data$salinity_conductivity),
                            max(meta_data$salinity_conductivity), 
                                by=0.1))

latitude.df <- meta_data %>%
  group_by(plant_species) %>%
  summarize(lat=median(lat),
            min_salinity=min(salinity_conductivity),
            max_salinity=max(salinity_conductivity)) %>%
  ungroup() %>%
  mutate(season="fall")

predictions.df <- left_join(newdata.df, latitude.df) %>%
  mutate(yi_pred = predict(gam.model, .)) %>%
  filter(between(salinity_conductivity, min_salinity, max_salinity))

ggplot() +
  geom_point(data = meta_data, aes(x = salinity_conductivity, y = sinh(yi), fill=plant_species, shape=plant_species), alpha = 0.5, size=2.5) +
  geom_line(data = predictions.df, aes(x = salinity_conductivity, y = sinh(yi_pred)), color = "black", size = 1) +
  scale_shape_manual(values = c("Cyperus malaccensis"=21, "Spartina alterniflora"=22,
                                "Distichlis Spicata" = 8, "Spartina patens" = 24, 
                                "Suaeda salsa" = 25, "Phragmites australis" = 23,
                                "Cladium jamaicense" = 13, "Juncae sp" = 20,
                                "Sesuvium portulacastrum" = 18, "Plantago maritima" = 10,
                                "Scirpus mariqueter" = 11, "Tamarix chinensis"=14,
                                "Salicornia sp."=12, "Suaeda sp." = 3))+
  scale_fill_manual(values = c("Cyperus malaccensis"="#025766", "Spartina patens" = "#5c3972", 
                               "Spartina alterniflora" = "#2FC45D", "Phragmites australis" ="#20bcff",
                               "Distichlis Spicata" = "#47724c", "Suaeda salsa" = "#c6b822",
                               "Cladium jamaicense" =  "#564e95", "Juncae sp" = "#86a291",
                               "Sesuvium portulacastrum" = "#90E669", "Plantago maritima" = "#47724c",
                               "Scirpus mariqueter" = "#60563f", "Tamarix chinensis"="#FF46A2",
                               "Salicornia sp."="black", "Suaeda sp." = "black")) +
  scale_y_continuous(transform = asinh_trans,
                     breaks=c(-1000,-100,-10,0,10,100,1000,10000))+
  geom_hline(yintercept=0, linetype="dashed", color="black")+
  facet_wrap(~ plant_species) +
  theme_classic() +
  theme(
    legend.position="none",
    axis.text.x = element_text(size = 11, colour = "black", family = "Arial"), 
    axis.text.y = element_text(size = 11, colour = "black", family = "Arial"),
    axis.title.x = element_text(size = 11, color = "black", family = "Arial"),
    axis.title.y = element_text(size = 1, color = "black", family = "Arial")) +
  labs(x = "Salinity", y =  expression(paste("CH"[4]*" flux (",mu,"mol m"^-2*" hr"^-1*")")))

ggsave("figures/Figure_S5.png",width=9,height=6)


